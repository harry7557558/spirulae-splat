# Code generation

Generated trees are marked `linguist-generated` / `linguist-vendored` in
`.gitattributes` and are **committed**. A fresh checkout therefore builds with
no Python at all — the dev build scripts run codegen when `python3` is
available and fall back to the committed files when it isn't.

Generated locations:

```
src/generated/       device-math headers from Slang
src/app/generated/   cli_config.h, viewer_html.h
src/backend/api/     backend module forwarders
src/instantiations/  kernel instantiation TUs
```

All generators run from the **repo root**:

```bash
python3 tools/codegen/generate_headers.py
python3 tools/codegen/generate_kernel_instantiation.py
python3 tools/codegen/generate_cli_config.py
python3 tools/codegen/generate_backend_api.py
# generate_vulkan_stubs.py takes arguments; see below
```

---

## `generate_headers.py` — `.cu` → `.cuh` declarations

Scans the `.cu` files listed in `HEADER_SOURCES` for functions preceded by the marker

```cpp
/*[AutoHeaderGeneratorExport]*/
void my_kernel_launcher(...) { ... }
```

and rewrites the declaration section of the matching `<Name>.cuh`.

**Invariants:**

1. Each `.cuh` has a splitter line:
   ```
   /* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */
   ```
   Everything **above** it is hand-written (includes, types, preamble) and
   preserved. Everything **below** is regenerated. Never edit below it.
2. `#if 0` blocks are stripped before scanning, so a disabled function is not
   exported.
3. **Sources are listed explicitly.** `HEADER_SOURCES` in the script maps each
   header stem to the ordered list of `.cu` files that declare into it. One
   header may be fed by many sources, and a source file's name is therefore
   free — name it after what it does. A listed file that does not exist raises
   `FileNotFoundError`, so a rename cannot silently drop declarations.
4. **A new kernel family needs an entry in `HEADER_SOURCES`.**
5. Declaration sections must stay CUDA-include-free — they must parse under
   `-DSSPLAT_BACKEND_VULKAN` with no CUDA toolkit present. Use
   `backend/api/BackendTypes.h` for vector types.

### Splitting a large `.cu` (the sanctioned workflow)

Because of invariant 3, splitting is nearly free:

```bash
# carve the sections out along the // ==== banners, naming each part after
# what it does -- not after the file it came from:
#   PixelWise.cu -> ImageConvert.cu, ImageColorOps.cu, DepthGeometry.cu,
#                   ImageDistort.cu, ImageWarp.cu, GtDepthNormalWarp.cu, Ppisp.cu
# put the shared preamble in <Family>Common.cuh and shared device helpers in a
# descriptively named header (BilinearSample.cuh), then:
#   - add the new files to HEADER_SOURCES['PixelWise']
python3 tools/codegen/generate_headers.py   # PixelWise.cuh's contents are unchanged
```

No build file changes: both build systems glob the kernel directories through
`cmake/sources.txt` (CMake with `CONFIGURE_DEPENDS`, so a new file is picked up
without a manual reconfigure). No caller changes: the header keeps its name and
its declaration set.

Naming rules when splitting:

- Name each part after its **function**, not its origin — `ImageWarp.cu`, not
  `PixelWise.Warp.cu`.
- `<Name>_kernel.cuh` is **reserved** for device bodies consumed by
  `generate_kernel_instantiation.py`. A plain shared-helper header gets a
  descriptive name instead.
- `<Family>Common.cuh` for the shared preamble (slang namespace includes,
  host-side view helpers).

---

## `generate_kernel_instantiation.py` — explicit template instantiations

Reads kernel declarations out of `src/kernels/**/*.cuh` and emits one small TU
per (primitive × camera model × SH degree × …) combination into
`src/instantiations/`
(currently 111 files). This keeps each nvcc edge small and parallelizable
rather than instantiating everything in one TU.

---

## `generate_cli_config.py` — Python dataclasses → native config

**The training config's single source of truth is the Python dataclasses** —
`TrainerConfig` and the nested `SpirulaeSplatDataParserConfig` /
`SpirulaeSplatDataManagerConfig` / `SpirulaeSplatModelConfig` /
`OptimizerConfig`, plus the tyro preset subclasses.

The script parses them with `ast` (no torch import — it runs on a fresh
checkout before `csrc.so` exists) and emits `src/app/generated/cli_config.h`:

- `struct SsplatConfig` — every field, flattened, with the Python defaults
  baked in;
- `SSPLAT_CONFIG_FIELDS(X)` — an X-macro over
  `(member, cli_key, group, choices, help)`, expanded by the CLI's generic
  parser and `--help` printer **and** by the GUI's "All Options" editor;
- `ssplat_apply_preset()` — one branch per preset, assigning exactly the
  fields that preset overrides.

Consequences: add a field to the Python dataclass and it appears in
`ssplat train --help`, in the GUI, and in the preset machinery after codegen.
Names are flattened (`--model.sh-degree` → `--sh-degree`); a collision across
groups must be resolved in the script's `RENAMES` table, and the script
**errors on any unlisted collision** so a new Python field cannot silently
shadow an existing flag.

Docstrings on the dataclass fields become CLI help text and GUI tooltips.

---

## `generate_backend_api.py` — backend module headers

Emits `src/backend/api/*.h` (`Projection.h`, `Rasterization.h`,
`PixelWise.h`, …) as **forwarders** that `#include` the relevant per-kernel
`.cuh` headers, grouped by subsystem.

Forwarding rather than copying declarations is deliberate: it keeps one
definition of every preamble type, so a TU may include both a module header
and an individual `.cuh` without ODR violations.

---

## `generate_vulkan_stubs.py` — link-probe for unported kernels

```bash
python3 tools/codegen/generate_vulkan_stubs.py <build-vulkan-dir> <output.cpp>
```

Link-probes `csrc_portable`'s objects against `libssplat_backend_vulkan.a`,
collects undefined C++ symbols, finds each function's declaration in the
per-kernel `.cuh` headers, and emits one **throwing** definition per function
(`backend/vulkan/kernels/TrainingStubs.gen.cpp`). This lets the portable
engine link today and lose stubs one module at a time as ports land.

Symbols with no header declaration (class methods, templates) are listed in
the generated file's header comment for manual stubbing in
`TrainingStubsManual.cpp`.

Rerun it whenever the engine gains a new kernel call, **or when a port lands**
— a ported symbol would otherwise be defined twice.

---

## Non-script generation done by CMake

- `viewer.html` → `app_generated/viewer_html.h` (byte array), so
  `ssplat train` and the GUI serve the viewer from a self-contained binary.
  Regenerated at configure time when the HTML changes. `Viewer.cpp` has a dev
  override so HTML edits don't require a rebuild.
- Slang → SPIR-V blobs → `vk_shaders_embedded.cpp`
  (`backend/vulkan/shaders/SpirvShaders.cmake` +
  `backend/vulkan/shaders/spirv_tool.cpp`). SPIR-V is **never committed**; one `slangc` edge per blob, with capability variants
  (`.atomicadd`, `.noint64`, `.int8`) compiled per entry point.
