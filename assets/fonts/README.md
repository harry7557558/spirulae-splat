# Fonts

Dear ImGui's built-in font is ASCII-only. Before this directory existed the
GUI could not draw a German umlaut, a French accent, a Turkish dotless i or a
single Cyrillic letter — which was already wrong, independently of
localization, for anyone whose dataset path was not pure ASCII.

Two faces, and they are handled differently for one reason: size.

| | file | size | how it ships |
|---|---|---|---|
| Latin / Cyrillic | `SpirulaUI-Regular.ttf` | 59 KB | **embedded** in the executable |
| CJK | `NotoSans{SC,TC,JP,KR}-Regular.otf` | 4–8 MB each | **fetched** at runtime, or bundled by a regional build |

## `SpirulaUI-Regular.ttf`

A subset of **Source Sans 3** (Adobe, SIL OFL 1.1 — `OFL-SourceSans.txt`)
covering exactly the ranges the twelve non-CJK languages in
`src/i18n/Languages.h` use: Basic Latin, Latin-1, Latin Extended-A, Cyrillic,
and the punctuation a translator will reach for. 431 KB upstream, 59 KB after
subsetting.

Source Sans 3 declares **"Source" as a Reserved Font Name**, so a Modified
Version may not carry it. Every name-table record was rewritten accordingly —
the copyright and licence records are quoted verbatim, as the OFL requires,
and everything else says *Spirula UI*. This is a licence obligation, not a
branding decision.

Rebuild it with:

```bash
pip install fonttools
python3 tools/make_ui_font.py           # downloads upstream, subsets, renames
python3 tools/make_ui_font.py --check   # verify the committed file
```

The output is byte-reproducible, so `--check` is meaningful. The committed
`.ttf` is a build artifact kept in the tree so that the build needs no network
— the same reasoning as the committed files under `src/generated/`.

## The CJK faces

Declared in `cjk_faces.txt`, which is the single source of truth: CMake
generates `app_generated/cjk_faces.h` from it and `src/app/gui/Fonts.cpp` reads
that. Nothing is committed here — they are downloaded, verified against the
SHA-256 in that file, and cached.

They are **Noto Sans CJK**, which is Source Han Sans under its other name, so
they are the same design as the Latin face above and a mixed line does not
visibly step. SIL OFL 1.1 with no reserved font name (`OFL-NotoSansCJK.txt`),
so they ship unmodified and under their own names.

There are four because **Han unification** gives the shared codepoints
different default glyph forms per region. A Japanese reader shown the
Simplified Chinese face gets kanji in Chinese forms: legible, and visibly
wrong. So the face follows the language, and `SS_FONT_CJK=sc` pairs with
`SS_DEFAULT_LANG=zh_hans`, not with an arbitrary one.

See `docs/i18n.md` for `SS_FONT_CJK`, the download prompt, and what a build
with `SS_FONT_CJK=none` does instead.
