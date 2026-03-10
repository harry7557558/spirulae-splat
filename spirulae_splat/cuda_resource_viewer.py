#!/usr/bin/env python3
"""
cuda_resource_viewer.py — Pretty-print CUDA kernel resource usage from compiled binaries.

Usage:
    python cuda_resource_viewer.py <path/to/binary.so>
    python cuda_resource_viewer.py <path/to/binary.so> --no-color
    python cuda_resource_viewer.py <path/to/binary.so> --sort-by reg|shared|const|name
"""

import sys
import re
import subprocess
import shutil
import argparse
from dataclasses import dataclass, field
from typing import Optional

# ---------------------------------------------------------------------------
# ANSI helpers
# ---------------------------------------------------------------------------

NO_COLOR = False  # toggled by --no-color


def _c(*codes: int) -> str:
    if NO_COLOR:
        return ""
    return f"\033[{';'.join(map(str, codes))}m"


RESET       = lambda: _c(0)
BOLD        = lambda: _c(1)
DIM         = lambda: _c(2)
ITALIC      = lambda: _c(3)
UNDERLINE   = lambda: _c(4)

# Foreground colours
BLACK   = lambda: _c(30)
RED     = lambda: _c(31)
GREEN   = lambda: _c(32)
YELLOW  = lambda: _c(33)
BLUE    = lambda: _c(34)
MAGENTA = lambda: _c(35)
CYAN    = lambda: _c(36)
WHITE   = lambda: _c(37)

# Bright variants
BRIGHT_RED     = lambda: _c(91)
BRIGHT_GREEN   = lambda: _c(92)
BRIGHT_YELLOW  = lambda: _c(93)
BRIGHT_BLUE    = lambda: _c(94)
BRIGHT_MAGENTA = lambda: _c(95)
BRIGHT_CYAN    = lambda: _c(96)
BRIGHT_WHITE   = lambda: _c(97)

# Background colours
BG_RED    = lambda: _c(41)
BG_YELLOW = lambda: _c(43)
BG_GREEN  = lambda: _c(42)
BG_BLUE   = lambda: _c(44)


def styled(text: str, *fns) -> str:
    prefix = "".join(f() for f in fns)
    if not prefix:
        return text
    return f"{prefix}{text}{RESET()}"


# ---------------------------------------------------------------------------
# Resource thresholds for colour-coding
# ---------------------------------------------------------------------------

REG_THRESHOLDS   = [(32, BRIGHT_GREEN), (64, BRIGHT_YELLOW), (96, BRIGHT_RED), (999, RED)]
SHARED_THRESHOLDS = [(0, DIM), (16384, BRIGHT_GREEN), (32768, BRIGHT_YELLOW), (999999, BRIGHT_RED)]


def color_value(value: int, thresholds: list) -> str:
    for limit, color_fn in thresholds:
        if value <= limit:
            return styled(str(value), color_fn)
    return str(value)


def reg_color(v: int) -> str:
    return color_value(v, REG_THRESHOLDS)


def shared_color(v: int) -> str:
    if v == 0:
        return styled("0", DIM)
    for limit, color_fn in SHARED_THRESHOLDS:
        if v <= limit:
            return styled(str(v), color_fn)
    return str(v)


def const_color(v: int) -> str:
    if v == 0:
        return styled("0", DIM)
    if v < 512:
        return styled(str(v), BRIGHT_GREEN)
    if v < 1024:
        return styled(str(v), BRIGHT_YELLOW)
    return styled(str(v), BRIGHT_RED)


def generic_color(v: int) -> str:
    if v == 0:
        return styled("0", DIM)
    return styled(str(v), BRIGHT_CYAN)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class KernelInfo:
    mangled_name: str
    demangled_name: Optional[str] = None
    reg: int = 0
    stack: int = 0
    shared: int = 0
    local: int = 0
    constants: dict = field(default_factory=dict)   # {idx: bytes}
    texture: int = 0
    surface: int = 0
    sampler: int = 0

    @property
    def display_name(self) -> str:
        return self.demangled_name or self.mangled_name


@dataclass
class FatbinSection:
    arch: str = ""
    code_version: str = ""
    host: str = ""
    compile_size: str = ""
    identifier: str = ""
    common_global: int = 0
    common_constants: dict = field(default_factory=dict)  # {idx: bytes}
    kernels: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

RES_PATTERN = re.compile(
    r"REG:(\d+)\s+STACK:(\d+)\s+SHARED:(\d+)\s+LOCAL:(\d+)"
    r"((?:\s+CONSTANT\[\d+\]:\d+)+)"
    r"\s+TEXTURE:(\d+)\s+SURFACE:(\d+)\s+SAMPLER:(\d+)"
)
CONST_PATTERN = re.compile(r"CONSTANT\[(\d+)\]:(\d+)")
COMMON_GLOBAL_PATTERN = re.compile(r"GLOBAL:(\d+)")
COMMON_CONST_PATTERN  = re.compile(r"CONSTANT\[(\d+)\]:(\d+)")


def parse_output(text: str) -> list[FatbinSection]:
    sections: list[FatbinSection] = []
    current: Optional[FatbinSection] = None
    in_resource_block = False
    pending_func: Optional[str] = None

    for line in text.splitlines():
        stripped = line.strip()

        # Detect new ELF fatbin section
        if stripped.startswith("Fatbin elf code:"):
            current = FatbinSection()
            sections.append(current)
            in_resource_block = False
            pending_func = None
            continue

        # Skip PTX sections entirely
        if stripped.startswith("Fatbin ptx code:"):
            current = None
            continue

        if current is None:
            continue

        # Header fields
        if stripped.startswith("arch ="):
            current.arch = stripped.split("=", 1)[1].strip()
        elif stripped.startswith("code version ="):
            current.code_version = stripped.split("=", 1)[1].strip()
        elif stripped.startswith("host ="):
            current.host = stripped.split("=", 1)[1].strip()
        elif stripped.startswith("compile_size ="):
            current.compile_size = stripped.split("=", 1)[1].strip()
        elif stripped.startswith("identifier ="):
            current.identifier = stripped.split("=", 1)[1].strip()
        elif stripped == "Resource usage:":
            in_resource_block = True
        elif not in_resource_block:
            continue
        elif stripped == "Common:":
            pass
        elif stripped.startswith("GLOBAL:"):
            m = COMMON_GLOBAL_PATTERN.search(stripped)
            if m:
                current.common_global = int(m.group(1))
            for cm in COMMON_CONST_PATTERN.finditer(stripped):
                current.common_constants[int(cm.group(1))] = int(cm.group(2))
        elif stripped.startswith("Function "):
            pending_func = stripped[len("Function "):]
            if pending_func.endswith(":"):
                pending_func = pending_func[:-1]
        elif pending_func and (stripped.startswith("REG:") or "REG:" in stripped):
            m = RES_PATTERN.search(stripped)
            if m:
                consts = {int(cm.group(1)): int(cm.group(2))
                          for cm in CONST_PATTERN.finditer(m.group(5))}
                k = KernelInfo(
                    mangled_name=pending_func,
                    reg=int(m.group(1)),
                    stack=int(m.group(2)),
                    shared=int(m.group(3)),
                    local=int(m.group(4)),
                    constants=consts,
                    texture=int(m.group(6)),
                    surface=int(m.group(7)),
                    sampler=int(m.group(8)),
                )
                current.kernels.append(k)
                pending_func = None

    return [s for s in sections if s.identifier or s.kernels]


# ---------------------------------------------------------------------------
# C++ demangling
# ---------------------------------------------------------------------------

def demangle_names(kernels: list[KernelInfo]) -> None:
    """Attempt to demangle mangled C++ names via c++filt."""
    if not kernels:
        return
    c_filt = shutil.which("c++filt")
    if not c_filt:
        return
    names = [k.mangled_name for k in kernels]
    try:
        result = subprocess.run(
            [c_filt] + names,
            capture_output=True, text=True, timeout=10
        )
        demangled = result.stdout.strip().splitlines()
        for k, d in zip(kernels, demangled):
            if d and d != k.mangled_name:
                k.demangled_name = d
    except Exception:
        pass


def demangle_all(sections: list[FatbinSection]) -> None:
    all_kernels = [k for s in sections for k in s.kernels]
    demangle_names(all_kernels)


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

BOX_H  = "─"
BOX_V  = "│"
BOX_TL = "╭"
BOX_TR = "╮"
BOX_BL = "╰"
BOX_BR = "╯"
BOX_LM = "├"
BOX_RM = "┤"

TERM_WIDTH = 100


def hr(char: str = BOX_H, width: int = TERM_WIDTH, color=None) -> str:
    line = char * width
    return styled(line, color) if color else line


def box_line(text: str, width: int = TERM_WIDTH) -> str:
    inner = width - 4
    return f" {styled(BOX_V, CYAN)}  {text:<{inner}}{styled(BOX_V, CYAN)}"


def section_header(title: str) -> str:
    pad = TERM_WIDTH - len(title) + 6
    left = pad // 2
    right = pad - left
    return (
        " " +
        styled(BOX_TL + BOX_H * (left + 1), CYAN) +
        styled(f" {title} ", BOLD, BRIGHT_WHITE) +
        styled(BOX_H * (right + 1) + BOX_TR, CYAN)
    )


def section_footer() -> str:
    return styled(" " + BOX_BL + BOX_H * (TERM_WIDTH - 3) + BOX_BR, CYAN)



def format_kernel_resources(k: KernelInfo) -> list[str]:
    """Return lines describing kernel resources (non-zero highlighted)."""
    items = []

    line_length = 0

    # REG is always shown
    items.append(f"{styled('REG', BOLD)}:{reg_color(k.reg)}")
    line_length += len('REG:') + len(str(k.reg))

    # STACK – show only if nonzero
    if k.stack:
        items.append(f"{styled('STACK', BOLD)}:{styled(str(k.stack), BRIGHT_RED)}")
        line_length += len('STACK:') + len(str(k.stack))

    # SHARED – always show (important for occupancy)
    if k.shared:
        items.append(f"{styled('SHARED', BOLD)}:{shared_color(k.shared)}")
        line_length += len('SHARED:') + len(str(k.shared))
    else:
        items.append(f"{styled('SHARED', DIM)}:{styled('0', DIM)}")
        line_length += len('SHARED:0')

    # LOCAL
    if k.local:
        items.append(f"{styled('LOCAL', BOLD)}:{BRIGHT_RED()}{k.local}{RESET()}")
        line_length += len('LOCAL:') + len(str(k.local))

    # CONSTANTS
    for idx in sorted(k.constants):
        v = k.constants[idx]
        label = styled(f"CONST[{idx}]", DIM)
        val   = const_color(v)
        items.append(f"{label}:{val}")
        line_length += len(f"CONST[{idx}]:") + len(str(val))

    # TEXTURE / SURFACE / SAMPLER (only if nonzero)
    for label, val in [("TEX", k.texture), ("SURF", k.surface), ("SAMP", k.sampler)]:
        if val:
            items.append(f"{styled(label, DIM)}:{styled(str(val), BRIGHT_MAGENTA)}")
            line_length += len(f"{label}:") + len(str(val))

    res_line = " ".join(items)
    line_length += len(items) - 1
    while line_length < 60:
        res_line += " "
        line_length += 1
    return res_line, line_length


def reg_bar(value: int, max_val: int = 255, width: int = 20) -> str:
    """Visual register usage bar."""
    if max_val == 0:
        return ""
    filled = round(value / max_val * width)
    filled = max(0, min(width, filled))
    ratio = value / max_val

    if ratio < 0.5:
        bar_color = BRIGHT_GREEN
    elif ratio < 0.75:
        bar_color = BRIGHT_YELLOW
    else:
        bar_color = BRIGHT_RED

    bar = styled("█" * filled, bar_color) + styled("░" * (width - filled), DIM)
    pct = f"{ratio*100:4.0f}%"
    return f"[{bar}] {styled(pct, DIM)}"


def truncate_name(name: str, max_len: int = TERM_WIDTH - 8) -> str:
    if len(name) <= max_len:
        return name + ' '*(max_len-len(name))
    if '>(' in name:
        name = name[:max(name.rfind('>(')+1, max_len-1)]
    else:
        name = name[:max_len - 1]
    return name + "…"


def print_section(section: FatbinSection, sort_key: Optional[str] = None, idx: int = 0) -> None:
    kernels = section.kernels
    if sort_key:
        if sort_key == "reg":
            kernels = sorted(kernels, key=lambda k: -k.reg)
        elif sort_key == "shared":
            kernels = sorted(kernels, key=lambda k: -k.shared)
        elif sort_key == "const":
            kernels = sorted(kernels, key=lambda k: -max(k.constants.values(), default=0))
        elif sort_key == "name":
            kernels = sorted(kernels, key=lambda k: k.display_name)

    # ---- Section header ----
    src_name = section.identifier.split("/")[-1] if "/" in section.identifier else section.identifier
    title = f"  {styled(src_name or f'Section {idx+1}', BOLD, BRIGHT_WHITE)}  "
    print()
    print(section_header(title.strip()))

    # ---- Metadata row ----
    meta_parts = [
        f"{styled('arch', DIM)}={styled(section.arch, BRIGHT_CYAN)}",
        f"{styled('ver', DIM)}={styled(section.code_version, BRIGHT_CYAN)}",
        f"{styled('host', DIM)}={styled(section.host, DIM)}",
        f"{styled('bits', DIM)}={styled(section.compile_size.replace('bit',''), DIM)}",
        f"{styled('kernels', DIM)}={styled(str(len(kernels)), BRIGHT_WHITE, BOLD)}",
    ]
    print(f" {styled('│', CYAN)}  " + "   ".join(meta_parts))

    # ---- Full identifier path ----
    if section.identifier:
        print(f" {styled('│', CYAN)}  {styled('src', DIM)}: {styled(section.identifier, DIM)}")

    # ---- Common resources ----
    common_parts = []
    if section.common_global:
        common_parts.append(f"{styled('GLOBAL', DIM)}:{generic_color(section.common_global)}")
    for idx2 in sorted(section.common_constants):
        v = section.common_constants[idx2]
        common_parts.append(f"{styled(f'CONST[{idx2}]', DIM)}:{const_color(v)}")
    if common_parts:
        print(f" {styled('│', CYAN)}  {styled('common', ITALIC, DIM)}: {' '.join(common_parts)}")

    # ---- Separator ----
    print(f" {styled(BOX_LM + BOX_H * (TERM_WIDTH - 3) + BOX_RM, CYAN)}")

    if not kernels:
        print(f" {styled('│', CYAN)}  {styled('(no kernels)', DIM)}")
    else:
        max_reg = max((k.reg for k in kernels), default=128)
        for i, k in enumerate(sorted(kernels, key=lambda k: k.display_name)):
            # Kernel name
            raw = k.display_name
            # Try to highlight template params dimly
            name_display = format_display_name(raw)
            print(f" {styled('│', CYAN)}  {styled(f'❯ {name_display}', BRIGHT_WHITE, BOLD)} {styled('│', CYAN)}")

            # Resources row
            res_line, line_length = format_kernel_resources(k)
            # bar = reg_bar(k.reg, max(max_reg, 32))
            bar = reg_bar(k.reg, 128)
            whitespace = ' ' * (TERM_WIDTH - line_length - 32)
            print(f" {styled('│', CYAN)}    {res_line}   {bar}  {whitespace} {styled('│', CYAN)}")

            # if i < len(kernels) - 1:
            #     print(f" {styled('│', CYAN)}  {styled('·' * (TERM_WIDTH - 6), DIM)} {styled('│', CYAN)}")

    print(section_footer())


def format_display_name(name: str) -> str:
    """Highlight template params and namespaces with dim styling."""
    name = truncate_name(name)

    # Color angle brackets for templates
    if NO_COLOR:
        return name

    result = ""
    depth = 0
    for ch in name:
        if ch == "<":
            depth += 1
            result += BRIGHT_MAGENTA() + ch + (DIM() if depth > 0 else "")
        elif ch == ">":
            depth -= 1
            result += (BRIGHT_MAGENTA() if depth == 0 else "") + ch + (RESET() + BRIGHT_WHITE() if depth == 0 else "")
        elif ch == ":" and depth == 0:
            result += DIM() + ch
        else:
            if depth > 0:
                result += DIM() + ch
            else:
                result += ch
    result += RESET()
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_cuobjdump(binary_path: str) -> str:
    tool = shutil.which("cuobjdump")
    if not tool:
        print(styled("✗ cuobjdump not found in PATH. Is CUDA toolkit installed?", BRIGHT_RED, BOLD),
              file=sys.stderr)
        sys.exit(1)

    try:
        result = subprocess.run(
            [tool, "--dump-resource-usage", binary_path],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            print(styled(f"✗ cuobjdump exited with code {result.returncode}:", BRIGHT_RED),
                  file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            sys.exit(result.returncode)
        return result.stdout
    except FileNotFoundError:
        print(styled("✗ cuobjdump not found.", BRIGHT_RED), file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print(styled("✗ cuobjdump timed out.", BRIGHT_RED), file=sys.stderr)
        sys.exit(1)


def print_summary(sections: list[FatbinSection]) -> None:
    total_kernels = sum(len(s.kernels) for s in sections)
    all_kernels = [k for s in sections for k in s.kernels]

    if not all_kernels:
        return

    max_reg_k    = max(all_kernels, key=lambda k: k.reg)
    max_shared_k = max(all_kernels, key=lambda k: k.shared)

    print()
    print(styled("  SUMMARY  ", BOLD, BRIGHT_WHITE, BG_BLUE) +
          styled(f"  {len(sections)} ELF section(s)   {total_kernels} kernel(s) total", DIM))
    print()

    # Per-architecture table
    archs: dict[str, list[KernelInfo]] = {}
    for s in sections:
        archs.setdefault(s.arch, []).extend(s.kernels)

    header = (
        styled(f"  {'ARCH':<10}", BOLD) +
        styled(f"{'KERNELS':>8}", BOLD) +
        styled(f"{'MAX REG':>10}", BOLD) +
        styled(f"{'AVG REG':>10}", BOLD) +
        styled(f"{'MAX SHARED':>12}", BOLD)
    )
    print(header)
    print(styled("  " + "─" * (len("ARCH") + 8 + 10 + 10 + 12 + 8), DIM))

    for arch, ks in sorted(archs.items()):
        avg_reg = sum(k.reg for k in ks) / len(ks) if ks else 0
        mx_reg  = max((k.reg for k in ks), default=0)
        mx_sh   = max((k.shared for k in ks), default=0)
        print(
            styled(f"  {arch:<10}", BRIGHT_CYAN) +
            f"{len(ks):>8}" +
            f"  {reg_color(mx_reg):>8}" + " " * 2 +
            f"  {styled(f'{avg_reg:.1f}', BRIGHT_WHITE):>8}" + " " * 2 +
            f"  {shared_color(mx_sh):>8}"
        )

    print()
    print(styled("  Highest REG:   ", DIM) +
          styled(f"{max_reg_k.reg}", BOLD, BRIGHT_RED) +
          "  " + styled(truncate_name(max_reg_k.display_name, TERM_WIDTH-40), DIM))
    if max_shared_k.shared > 0:
        print(styled("  Highest SHARED:", DIM) +
              styled(f" {max_shared_k.shared}", BOLD, BRIGHT_YELLOW) +
              "  " + styled(truncate_name(max_shared_k.display_name, TERM_WIDTH-40), DIM))
    print()


def main() -> None:
    global NO_COLOR

    parser = argparse.ArgumentParser(
        description="Pretty-print CUDA kernel resource usage from a compiled binary.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("binary", help="Path to CUDA binary (.so, .dll, …)")
    parser.add_argument("--no-color", action="store_true", help="Disable ANSI colour output")
    parser.add_argument(
        "--sort-by",
        choices=["reg", "shared", "const", "name"],
        default=None,
        help="Sort kernels within each section",
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        help="Print raw cuobjdump output and exit (for debugging)",
    )
    args = parser.parse_args()

    NO_COLOR = args.no_color or not sys.stdout.isatty()

    raw = run_cuobjdump(args.binary)

    if args.raw:
        print(raw)
        return

    sections = parse_output(raw)
    demangle_all(sections)

    if not sections:
        print(styled("No ELF fatbin sections with resource usage found.", BRIGHT_YELLOW))
        sys.exit(0)

    # Banner
    print()
    print(styled("━" * TERM_WIDTH, CYAN))
    print(styled("  CUDA KERNEL RESOURCE VIEWER", BOLD, BRIGHT_CYAN) +
          styled(f"  ·  {args.binary}", DIM))
    print(styled("━" * TERM_WIDTH, CYAN))

    print_summary(sections)

    for i, section in enumerate(sections):
        print_section(section, sort_key=args.sort_by, idx=i)

    print()


if __name__ == "__main__":
    main()
