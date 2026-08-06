#!/usr/bin/env python3
"""Rebuild the five font files the GUI embeds.

    assets/fonts/SpirulaUI-Regular.ttf    Latin + Cyrillic, from Source Sans 3
    assets/fonts/SpirulaCJK-JP.otf   \
    assets/fonts/SpirulaCJK-SC.otf    |   from Noto Sans CJK, subset to the
    assets/fonts/SpirulaCJK-TC.otf    |   characters this program's own
    assets/fonts/SpirulaCJK-KR.otf   /    translations actually use

Together they are ~490 KB and they are what makes a default build render all
thirteen languages in SS_LANGUAGES with nothing to download. That is worth
being precise about, because the obvious alternative is not:

  * A full Noto Sans CJK face is 4-8 MB, and there are four of them, because
    Han unification gives the shared codepoints different default glyph forms
    per region. Embedding all four is 23 MB of executable to render a menu.
  * Downloading one on demand -- what this GUI did before -- means the first
    thing a Japanese user sees after picking their language is a wall of
    boxes and a download button. The language picker is exactly where a user
    who cannot read the current UI language has arrived, so it is exactly the
    wrong place to require reading.

Subsetting to the ~600 characters per region the UI actually writes costs
~110 KB per face, so all four fit. The full faces are still fetched on demand
(src/app/gui/Fonts.h) -- not for the UI, but for dataset paths and file names,
which are user data and can hold any character at all.

THE SUBSETS ARE DERIVED FROM THE CATALOGS. Edit a translation and they are
stale, which would show up as one boxed character in the middle of a sentence.
tools/check_font_coverage.py is the guard: it runs on every build, needs no
network and no fontTools, and fails when a catalog uses a character no
committed font has.

Two more things this script exists to keep honest:

  * A subset is a "Modified Version" under the SIL Open Font License. Source
    Sans 3 reserves the name "Source", so every name-table record carrying it
    has to be rewritten -- not just the filename. Noto reserves nothing, but
    a 600-glyph file called "Noto Sans JP" would still be a lie, so it is
    renamed too. Copyright (nameID 0) and licence (13) are quoted verbatim,
    as the licence requires.
  * The committed files are build artifacts. This is their generator, in the
    same spirit as tools/codegen/ -- the output is committed so the build
    needs no network, and this script is how you check it is what it claims.

Usage:
    pip install fonttools
    python3 tools/make_ui_font.py            # rebuild all five
    python3 tools/make_ui_font.py --check    # verify the committed files
    python3 tools/make_ui_font.py --latin    # just the Latin face
"""

import argparse
import hashlib
import io
import pathlib
import re
import subprocess
import sys
import tempfile
import urllib.request
import zipfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
FONT_DIR = ROOT / "assets" / "fonts"
CATALOG_DIR = ROOT / "src" / "i18n" / "catalog"
LANGUAGES_H = ROOT / "src" / "i18n" / "Languages.h"
CJK_SPEC = FONT_DIR / "cjk_faces.txt"

SOURCE_SANS_URL = (
    "https://github.com/adobe-fonts/source-sans/releases/download/"
    "3.052R/TTF-source-sans-3.052R.zip"
)
SOURCE_SANS_MEMBER = "TTF/SourceSans3-Regular.ttf"

# What the twelve non-CJK locales in SS_LANGUAGES actually need. Ranges, not
# a character list, because a translator will reach for a dash or a quote mark
# nobody enumerated and tofu in the middle of a sentence is worse than 3 KB.
UNICODES = ",".join([
    "U+0020-007E",     # Basic Latin
    "U+00A0-00FF",     # Latin-1: de/fr/es/pt/it/nl accents
    "U+0100-017F",     # Latin Extended-A: Turkish g-breve, s-cedilla, dotted I
    "U+0192",          # florin, occasionally used as a function sign
    "U+01FA-01FF",
    "U+02C6-02DD",     # spacing modifiers (circumflex, caron, ...)
    "U+0394,U+03A9,U+03BC,U+03C0",   # delta/omega/mu/pi -- they appear in help text
    "U+0400-045F",     # Cyrillic: Russian
    "U+0490-0491",     # Ukrainian ghe -- free, and stops one obvious hole
    "U+2000-206F",     # dashes, curly quotes, ellipsis, bullet
    "U+20AC",          # euro
    "U+2122",          # trademark
    "U+2190-2193",     # arrows -- the UI writes "-> outputs/scene"
    "U+2202,U+2206,U+2211,U+2212,U+221A,U+221E,U+2248,U+2260,U+2264,U+2265",
    "U+25A0,U+25CF",   # square/circle bullets
    "U+FB01-FB02",     # fi/fl ligatures
    "U+FFFD",          # replacement character -- what a bad byte should look like
])

LAYOUT_FEATURES = "kern,liga,ccmp,locl,mark,mkmk"

# The catalog tag macro (src/i18n/BeginCatalog.h) -> the regional face whose
# glyph forms that language is read in.
TAG_REGION = {"JA": "jp", "ZH_HANS": "sc", "ZH_HANT": "tc", "KO": "kr"}

# Characters every regional subset carries whether or not that language's
# translations use them.
#
#   - the native name of every language in SS_LANGUAGES, because the language
#     picker draws all thirteen at once and the whole point of this exercise
#     is that a user can read their own language's name in it
#   - SS_LANG_MENU_ICON, the language menu's label
#   - CJK punctuation, because the alternative to two dozen free glyphs is
#     regenerating four fonts because a translator wrote a full-width comma
PUNCTUATION = "、。，．：；？！「」『』（）〔〕【】・…‥ー～〜／＼％＋－＝＜＞　"


# ---------------------------------------------------------------------------
# Scanning the catalogs
# ---------------------------------------------------------------------------

# TAG("..." "...") -- adjacent string literals, possibly across lines.
_LITERALS = r'(?:"(?:[^"\\]|\\.)*"\s*)+'
_ONE_LITERAL = re.compile(r'"((?:[^"\\]|\\.)*)"')
_ESCAPES = {"n": "\n", "t": "\t", "r": "\r", "0": "\0",
            '"': '"', "\\": "\\", "'": "'"}


def _unescape(blob: str) -> str:
    """Concatenated C string literals -> the text they denote."""
    out = []
    for m in _ONE_LITERAL.finditer(blob):
        s, i = m.group(1), 0
        while i < len(s):
            if s[i] == "\\" and i + 1 < len(s):
                out.append(_ESCAPES.get(s[i + 1], s[i + 1]))
                i += 2
            else:
                out.append(s[i])
                i += 1
    return "".join(out)


def _languages_h() -> list:
    """SS_LANGUAGES -> [(id, code, native, english)], the one source of truth."""
    text = LANGUAGES_H.read_text(encoding="utf-8")
    rows = re.findall(
        r'X\(\s*(\w+)\s*,\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,\s*"([^"]*)"\s*\)', text)
    if len(rows) < 2:
        raise SystemExit(f"could not read SS_LANGUAGES out of {LANGUAGES_H}")
    return rows


def native_names() -> str:
    """The `native` column of SS_LANGUAGES -- what the picker draws."""
    return "".join(r[2] for r in _languages_h())


def always() -> set:
    """What every regional subset carries regardless of its translations."""
    text = LANGUAGES_H.read_text(encoding="utf-8")
    m = re.search(r'#define\s+SS_LANG_MENU_ICON\s+"([^"]*)"', text)
    if not m:
        raise SystemExit(f"could not read SS_LANG_MENU_ICON out of {LANGUAGES_H}")
    return set(PUNCTUATION) | set(m.group(1)) | set(native_names())


def _tag_hits(text: str, tag: str):
    """Every TAG("...") in `text`, unescaped."""
    for m in re.finditer(tag + r'\(\s*(' + _LITERALS + r')\)', text):
        yield _unescape(m.group(1))


def _catalog_text() -> str:
    return "\n".join(p.read_text(encoding="utf-8")
                     for p in sorted(CATALOG_DIR.glob("*.h")))


def scan_catalogs() -> dict:
    """region id -> the characters that region's translations use.

    Shared with tools/check_font_coverage.py, which imports it.
    """
    per = {r: set(always()) for r in TAG_REGION.values()}
    text = _catalog_text()
    for tag, region in TAG_REGION.items():
        for s in _tag_hits(text, tag):
            per[region] |= set(s)
    # ASCII is the Latin face's job; a regional face carrying its own Latin
    # would win the merge for it and quietly restyle every English string.
    return {r: {c for c in cs if ord(c) > 0x7F} for r, cs in per.items()}


def scan_all() -> set:
    """Every character the UI can draw from its own catalogs, all languages.

    The union the committed fonts have to cover between them. Tags are derived
    from SS_LANGUAGES (`zh_hans` -> `ZH_HANS`), so a new language is scanned
    the moment it exists rather than when someone remembers this file.
    """
    chars = always()
    text = _catalog_text()
    for row in _languages_h():
        for s in _tag_hits(text, row[0].upper()):
            chars |= set(s)
    return chars


def read_cjk_spec() -> list:
    """assets/fonts/cjk_faces.txt -> [{id, file, url, sha256, bytes, label}]."""
    faces = []
    for line in CJK_SPEC.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split(" | ")
        if len(f) != 6:
            raise SystemExit(f"{CJK_SPEC}: expected 6 fields, got {len(f)}:\n  {line}")
        faces.append(dict(zip(("id", "file", "url", "sha256", "bytes", "label"), f)))
    return faces


# ---------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------

def fetch(url: str, expect_sha: str = "") -> bytes:
    print(f"downloading {url}")
    with urllib.request.urlopen(url, timeout=300) as r:
        blob = r.read()
    if expect_sha:
        got = hashlib.sha256(blob).hexdigest()
        if got != expect_sha:
            raise SystemExit(f"{url}\n  expected sha256 {expect_sha}\n  got      {got}")
    return blob


def rename(font, *, family: str, ps_name: str, version: str,
           description: str, forbidden: str = "") -> None:
    """Give the subset its own identity, as the OFL requires of a modification."""
    name = font["name"]
    replacements = {
        1: family,
        2: "Regular",
        3: f"{ps_name};{version}",
        4: family,
        5: version,
        6: ps_name,
        10: description,
        16: family,
        17: "Regular",
    }
    for record in list(name.names):
        # 7 is the trademark notice, which belongs to the upstream foundry and
        # must not be reused; 8/9/11/12/19-22/25 are vendor and sample strings
        # that would all be claims about someone else.
        if record.nameID in (7, 8, 9, 11, 12, 19, 20, 21, 22, 25):
            name.names.remove(record)
            continue
        new = replacements.get(record.nameID)
        if new is not None:
            record.string = new
    if forbidden:
        for record in name.names:
            if record.nameID in (0, 13):
                continue      # copyright + licence: quoted verbatim, as required
            if forbidden in str(record):
                raise SystemExit(
                    f"nameID {record.nameID} still says '{forbidden}': {record}")


def _subset(src: pathlib.Path, dst: pathlib.Path, args: list) -> None:
    subprocess.run(
        [sys.executable, "-m", "fontTools.subset", str(src),
         "--no-hinting", "--desubroutinize", "--drop-tables+=DSIG",
         f"--output-file={dst}"] + args,
        check=True)


def build_latin() -> bytes:
    from fontTools.ttLib import TTFont

    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        src, dst = tmp / "in.ttf", tmp / "out.ttf"
        with zipfile.ZipFile(io.BytesIO(fetch(SOURCE_SANS_URL))) as z:
            src.write_bytes(z.read(SOURCE_SANS_MEMBER))
        _subset(src, dst, [f"--unicodes={UNICODES}",
                           f"--layout-features={LAYOUT_FEATURES}"])
        # recalcTimestamp=False, or head.modified is stamped with the wall
        # clock on save and the committed file differs on every run -- which
        # would make --check useless.
        font = TTFont(dst, recalcTimestamp=False)
        rename(font, family="Spirula UI", ps_name="SpirulaUI-Regular",
               version="Version 3.052;Spirula subset",
               description=(
                   "Subset of Source Sans 3 (Adobe, SIL OFL 1.1) covering the "
                   "Latin and Cyrillic ranges Spirula Studio's user interface "
                   "uses. Renamed as the OFL requires of a Modified Version: "
                   "'Source' is a Reserved Font Name."),
               forbidden="Source")
        out = tmp / "final.ttf"
        font.save(out)
        return out.read_bytes()


def build_cjk(face: dict, chars: set) -> bytes:
    from fontTools.ttLib import TTFont

    region = face["id"].upper()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        src, dst = tmp / "in.otf", tmp / "out.otf"
        src.write_bytes(fetch(face["url"], face["sha256"]))

        # Which of the wanted characters this regional face even has. A
        # simplified-only character has no Japanese form to draw, so it is
        # absent here and comes from the SC subset at merge time -- see
        # Fonts.cpp. Reported rather than fatal: it is a translation smell,
        # not a build error.
        have = set(TTFont(src).getBestCmap())
        missing = {c for c in chars if ord(c) not in have}
        if missing:
            print(f"  note: {face['file']} has no glyph for "
                  f"{''.join(sorted(missing))} -- another face will supply it")

        # No layout features: CJK UI text needs no kerning or shaping, and
        # GSUB/GPOS for a CJK face is most of its weight.
        codes = ",".join(f"U+{ord(c):04X}" for c in sorted(chars - missing))
        _subset(src, dst, [f"--unicodes={codes}", "--layout-features="])
        font = TTFont(dst, recalcTimestamp=False)
        rename(font, family=f"Spirula CJK {region}",
               ps_name=f"SpirulaCJK-{region}",
               version="Version 2.004;Spirula subset",
               description=(
                   f"Subset of Noto Sans CJK {region} (Google, SIL OFL 1.1) "
                   f"covering only the characters Spirula Studio's own "
                   f"{face['label']} translations use. Renamed because a file "
                   f"this small is not the font it was cut from."))
        out = tmp / "final.otf"
        font.save(out)
        return out.read_bytes()


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def latin_coverage(blob: bytes) -> None:
    from fontTools.ttLib import TTFont

    cmap = TTFont(io.BytesIO(blob)).getBestCmap()
    samples = {
        "German": "äöüßÄÖÜ",
        "French": "éèêëçàùîïôœÉÀ«»",
        "Spanish": "ñáíóúü¿¡",
        "Portuguese": "ãõâêçÁÍ",
        "Italian": "àèéìòù",
        "Dutch": "ëïĳ",
        "Turkish": "ğĞşŞıİçÇöÖüÜ",
        "Russian": ("АБВГДЕЁЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯ"
                    "абвгдеёжзийклмнопрстуфхцчшщъыьэюя"),
        "punctuation": "–—‘’“”…•·×°→←",
    }
    bad = False
    for label, chars in samples.items():
        missing = "".join(c for c in chars if ord(c) not in cmap)
        print(f"  {label:12s} {'OK' if not missing else 'MISSING ' + missing}")
        bad = bad or bool(missing)
    if bad:
        raise SystemExit("coverage check failed")


def report(path: pathlib.Path, blob: bytes, glyphs: int = 0) -> None:
    print(f"  {path.name:24s} {len(blob) / 1024:6.0f} KB"
          + (f"  {glyphs} characters" if glyphs else ""))


# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="rebuild and compare against the committed files")
    ap.add_argument("--latin", action="store_true", help="only the Latin face")
    ap.add_argument("--cjk", action="store_true", help="only the CJK subsets")
    args = ap.parse_args()
    do_latin = args.latin or not args.cjk
    do_cjk = args.cjk or not args.latin

    built = {}
    if do_latin:
        blob = build_latin()
        latin_coverage(blob)
        built[FONT_DIR / "SpirulaUI-Regular.ttf"] = blob
        report(FONT_DIR / "SpirulaUI-Regular.ttf", blob)

    if do_cjk:
        per_region = scan_catalogs()
        for face in read_cjk_spec():
            chars = per_region[face["id"]]
            blob = build_cjk(face, chars)
            path = FONT_DIR / f"SpirulaCJK-{face['id'].upper()}.otf"
            built[path] = blob
            report(path, blob, len(chars))

    total = sum(len(b) for b in built.values())
    print(f"  {'total':24s} {total / 1024:6.0f} KB")

    if args.check:
        bad = False
        for path, blob in built.items():
            want = hashlib.sha256(blob).hexdigest()
            have = (hashlib.sha256(path.read_bytes()).hexdigest()
                    if path.exists() else "(missing)")
            if have != want:
                print(f"FAIL {path.relative_to(ROOT)}\n"
                      f"  committed {have}\n  rebuilt   {want}")
                bad = True
        if bad:
            raise SystemExit("committed fonts differ from a fresh build; "
                             "re-run without --check")
        print(f"OK: {len(built)} committed font(s) match a fresh build")
        return

    FONT_DIR.mkdir(parents=True, exist_ok=True)
    for path, blob in built.items():
        path.write_bytes(blob)
        print(f"wrote {path.relative_to(ROOT)}  "
              f"sha256 {hashlib.sha256(blob).hexdigest()}")


if __name__ == "__main__":
    main()
