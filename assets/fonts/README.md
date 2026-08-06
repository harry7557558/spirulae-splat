# Fonts

Dear ImGui's built-in font is ASCII-only. Before this directory existed the
GUI could not draw a German umlaut, a French accent, a Turkish dotless i or a
single Cyrillic letter — which was already wrong, independently of
localization, for anyone whose dataset path was not pure ASCII.

Five faces are committed here and **all five are embedded** in the executable,
507 KB together. Nothing has to be downloaded for the interface to render in
any of the thirteen languages in `src/i18n/Languages.h`.

| | file | size | from |
|---|---|---|---|
| Latin / Cyrillic | `SpirulaUI-Regular.ttf` | 58 KB | Source Sans 3 |
| CJK, per region | `SpirulaCJK-{JP,SC,TC,KR}.otf` | 449 KB | Noto Sans CJK |

The full CJK faces — 4–8 MB each — are still **fetched at runtime** or bundled
by a regional build, but for a much smaller job: see *The full faces* below.

## `SpirulaUI-Regular.ttf`

A subset of **Source Sans 3** (Adobe, SIL OFL 1.1 — `OFL-SourceSans.txt`)
covering exactly the ranges the twelve non-CJK languages use: Basic Latin,
Latin-1, Latin Extended-A, Cyrillic, and the punctuation a translator will
reach for. 431 KB upstream, 58 KB after subsetting.

Source Sans 3 declares **"Source" as a Reserved Font Name**, so a Modified
Version may not carry it. Every name-table record was rewritten accordingly —
the copyright and licence records are quoted verbatim, as the OFL requires,
and everything else says *Spirula UI*. This is a licence obligation, not a
branding decision.

## `SpirulaCJK-{JP,SC,TC,KR}.otf`

Each is **Noto Sans CJK** (Google, SIL OFL 1.1 — `OFL-NotoSansCJK.txt`) cut
down to the ~600 characters that region's own translations use, plus the
native name of every language and a little CJK punctuation. ~110 KB each.
Noto Sans CJK is Source Han Sans under its other name, so it is the same
design as the Latin face and a mixed line does not visibly step.

**There are four because of Han unification**: the shared codepoints have
different default glyph forms per region, so a Japanese reader shown the
Simplified Chinese face gets kanji in Chinese forms — legible, and visibly
wrong. Two alternatives were measured and rejected:

| approach | size | cost |
|---|---|---|
| four regional subsets *(what ships)* | 449 KB | — |
| one pan-CJK subset | 287 KB | 262 shared characters in the wrong regional form for three of the four languages |
| a shared base + three deltas | 358 KB | of the characters used by more than one region, only 267 have identical outlines and 262 genuinely differ — so the deltas are most of the weight anyway |

Noto reserves no font name, so a subset could legally keep it; these are
renamed anyway, because a 600-glyph file called "Noto Sans JP" is a lie.

## Rebuilding

```bash
pip install fonttools
python3 tools/make_ui_font.py           # all five: download, subset, rename
python3 tools/make_ui_font.py --check   # verify the committed files
python3 tools/check_font_coverage.py    # cheap: did a translation outgrow them?
```

The output is byte-reproducible, so `--check` is meaningful. The committed
files are build artifacts kept in the tree so that the build needs no network
— the same reasoning as the committed files under `src/generated/`.

**The CJK subsets are derived from the catalogs.** Edit a translation and they
go stale, and the symptom is one hollow box in the middle of an otherwise fine
sentence. `tools/check_font_coverage.py` is the guard and runs on every
`build_develop.bash`; it needs neither the network nor fontTools, because it
parses `cmap` by hand.

## The full faces

Declared in `cjk_faces.txt`, which is the single source of truth: CMake
generates `app_generated/cjk_faces.h` from it and `src/app/gui/Fonts.cpp` reads
that. Nothing of them is committed here — they are downloaded, verified against
the SHA-256 in that file, and cached.

They exist for the text this program did **not** write: dataset paths, file
names and typed mask prompts are user data and can hold any character at all,
which no subset of our own strings can anticipate. Until one is fetched, a
folder called `C:\写真\` may render as boxes; the interface itself never does.

See `docs/i18n.md` for `SS_FONT_CJK`, the download offer, and what a build with
`SS_FONT_CJK=none` does instead.
