#!/usr/bin/env python3
"""Compile Vulkan compute entry points (slang/vulkan/*.slang) to SPIR-V.

The SPIR-V counterpart of slang/build.bash, but the blobs are NOT committed:
CMake (SSPLAT_BACKEND=vulkan) runs this at build time into the build tree
and embeds the result via spirulae_splat/embed_spirv.py. CMake locates
slangc (or downloads the pinned release) and passes it via --slangc.

Entry points are discovered by scanning for [shader("compute")] markers and
line-initial DEF_*(name, ...) macro instantiations. Each entry compiles to
its own <stem>.<entry>.spv. Entries whose include closure calls
atomic_add_f32 (see vulkan/atomic_float.slang) additionally compile a
<stem>.<entry>.atomicadd.spv variant with -DSSPLAT_NATIVE_F32_ATOMIC_ADD
(native OpAtomicFAddEXT instead of the CAS-loop emulation); the pipeline
layer (VulkanPipelines.cpp) picks that blob on devices with
shaderBufferFloat32AtomicAdd. A checksum cache (<out-dir>/checksums.json)
keyed per entry over the source, its transitive #include closure, the
compile arguments, and this script skips up-to-date blobs.

Usage: python3 slang/build_spirv.py --out-dir DIR [--force] [--slangc PATH]
(run from the `cuda` folder)
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

SRC_DIR = "slang/vulkan"
INCLUDE_DIRS = ["slang", "slang/vulkan"]

COMPILE_ARGS = [
    "-target", "spirv",
    "-stage", "compute",
    "-O2",
    # "-ignore-capabilities",
    # "-line-directive-mode", "none",
    "-g2",
]

# Native float-atomic variant (must match kNativeF32AtomicSuffix in
# csrc/backend/vulkan/VulkanPipelines.cpp).
ATOMIC_HEADER = os.path.normpath(os.path.join(SRC_DIR, "atomic_float.slang"))
ATOMIC_SUFFIX = ".atomicadd"
ATOMIC_DEFINE = "-DSSPLAT_NATIVE_F32_ATOMIC_ADD"
ATOMIC_CALL_RE = re.compile(r'\batomic_add_f32\s*\(')

ENTRY_RE = re.compile(
    r'\[shader\("compute"\)\]\s*(?:(?:\[[^\]]*\]|//[^\n]*)\s*)*void\s+(\w+)\s*\(')
# Macro-instantiated entry points: a line-initial `DEF_FOO(entry_name, ...)`
# where DEF_FOO is a local macro wrapping [shader("compute")].
MACRO_ENTRY_RE = re.compile(r'^DEF_\w+\(\s*(\w+)\s*[,)]', re.MULTILINE)
INCLUDE_RE = re.compile(r'^\s*#include\s+"([^"]+)"', re.MULTILINE)


def find_entries(path):
    src = open(path).read()
    # Both regexes can capture the literal macro parameter "NAME" from a
    # #define body; instantiations use real identifiers.
    seen = set()
    out = []
    for name in ENTRY_RE.findall(src) + MACRO_ENTRY_RE.findall(src):
        if name != "NAME" and name not in seen:
            seen.add(name)
            out.append(name)
    return out


def include_closure(path, seen=None):
    """Transitive #include closure, resolved against INCLUDE_DIRS and the
    including file's directory."""
    if seen is None:
        seen = []
    if path in seen or not os.path.exists(path):
        return seen
    seen.append(path)
    src = open(path).read()
    for inc in INCLUDE_RE.findall(src):
        for base in [os.path.dirname(path)] + INCLUDE_DIRS:
            cand = os.path.normpath(os.path.join(base, inc))
            if os.path.exists(cand):
                include_closure(cand, seen)
                break
    return seen


def entry_checksum(source, entry, args):
    h = hashlib.blake2s()
    for dep in sorted(include_closure(source)):
        h.update(dep.encode())
        h.update(open(dep, "rb").read())
    h.update(entry.encode())
    h.update(" ".join(args).encode())
    h.update(open(__file__, "rb").read())
    return h.hexdigest()


def uses_float_atomic_add(source):
    """True when the entry's include closure calls atomic_add_f32 (the
    definition in atomic_float.slang itself does not count)."""
    closure = include_closure(source)
    if ATOMIC_HEADER not in (os.path.normpath(p) for p in closure):
        return False
    return any(ATOMIC_CALL_RE.search(open(dep).read())
               for dep in closure if os.path.normpath(dep) != ATOMIC_HEADER)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--slangc", default=None)
    opts = ap.parse_args()
    out_dir = opts.out_dir
    cache_file = os.path.join(out_dir, "checksums.json")

    slangc = opts.slangc or shutil.which("slangc")
    if not slangc:
        sys.exit("slangc not found in PATH (tested with slang-2026.2.1); "
                 "pass --slangc")

    os.makedirs(out_dir, exist_ok=True)
    cache = {}
    if os.path.exists(cache_file) and not opts.force:
        try:
            cache = json.load(open(cache_file))
        except (json.JSONDecodeError, OSError):
            cache = {}

    sources = [os.path.join(SRC_DIR, f) for f in sorted(os.listdir(SRC_DIR))
               if f.endswith(".slang")]
    # Files pulled in as headers by another source; entry-less .slang files
    # are only suspicious when nothing includes them.
    included = set()
    for source in sources:
        included.update(os.path.normpath(p) for p in include_closure(source)
                        if os.path.normpath(p) != os.path.normpath(source))

    jobs = []
    expected = set()  # every (stem, entry, variant) key live in the sources
    for source in sources:
        fname = os.path.basename(source)
        stem = fname[:-len(".slang")]
        entries = find_entries(source)
        if not entries and os.path.normpath(source) not in included:
            print(f"warning: no compute entry points in {source}")
        variants = [("", [])]
        if uses_float_atomic_add(source):
            variants.append((ATOMIC_SUFFIX, [ATOMIC_DEFINE]))
        for entry in entries:
            for suffix, extra_args in variants:
                key = f"{stem}.{entry}{suffix}"
                expected.add(key)
                out = os.path.join(out_dir, key + ".spv")
                args = COMPILE_ARGS + extra_args
                checksum = entry_checksum(source, entry, args)
                if cache.get(key) == checksum and os.path.exists(out):
                    print(f"  up-to-date: {out}")
                    continue
                jobs.append((source, entry, out, key, checksum, args))

    def compile_one(job):
        source, entry, out, key, checksum, args = job
        cmd = [slangc, source, "-entry", entry, "-o", out]
        cmd += [f"-I{d}" for d in INCLUDE_DIRS]
        cmd += args
        proc = subprocess.run(cmd, capture_output=True, text=True)
        return key, checksum, out, proc

    failed = False
    with ThreadPoolExecutor() as pool:
        for key, checksum, out, proc in pool.map(compile_one, jobs):
            if proc.returncode != 0:
                print(f"  FAILED: {out}\n{proc.stdout}{proc.stderr}",
                      file=sys.stderr)
                failed = True
                cache.pop(key, None)
            else:
                if proc.stderr.strip():
                    print(f"  (slangc warnings for {out})\n{proc.stderr}")
                print(f"  built: {out}")
                cache[key] = checksum

    # Drop cache entries and blobs whose entry (or variant) no longer exists.
    cache = {k: v for k, v in cache.items() if k in expected}
    live = {k + ".spv" for k in expected}
    for f in os.listdir(out_dir):
        if f.endswith(".spv") and f not in live:
            os.remove(os.path.join(out_dir, f))
            print(f"  removed stale: {f}")

    json.dump(cache, open(cache_file, "w"), indent=1, sort_keys=True)
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
