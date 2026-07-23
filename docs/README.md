# Documentation index

Start at [`../AGENTS.md`](../AGENTS.md) for orientation. This directory holds
the detail.

| document | what it covers |
|---|---|
| [architecture.md](architecture.md) | layers, data flow, engine state, where each responsibility lives |
| [build.md](build.md) | the build matrix, every option, per-platform notes |
| [backends.md](backends.md) | the CUDA/Vulkan seam, coverage, pointers to the authoritative backend docs |
| [codegen.md](codegen.md) | the five generators, what they own, and the invariants |
| [datasets.md](datasets.md) | supported dataset layouts and how they are parsed |
| [testing.md](testing.md) | native parity tests, the three Python/C++ parity gates, the reference-dump workflow |
| [restructure-proposal.md](restructure-proposal.md) | the in-progress plan for reorganizing the tree |
| [notes/pose-normalization.md](notes/pose-normalization.md) | orientation/centering: what the native parser implements, and the kept Python reference for what it doesn't |
| [notes/](notes/) | design notes for individual subsystems |

Authoritative documents that live next to their code rather than here:

- `src/backend/README.md` — the backend seam design.
- `src/backend/vulkan/README.md` — the Vulkan
  backend in depth: device baseline, capability variants, memory model,
  atomics, Slang notes, kernel coverage. **The single most detailed document
  in the repo**; read it before touching Vulkan code.
- `src/app/README.md` — working notes on the
  standalone CLI and the native GUI.
- `viewer/README.md` — the standalone WebGL2/WASM viewer.
- `scripts/README.md` — dataset preprocessing tools.
