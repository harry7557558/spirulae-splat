#!/usr/bin/env python3
"""Rebuild assets/fonts/SpirulaUI-Regular.ttf from Source Sans 3.

The GUI embeds one Latin/Cyrillic face so that a default build renders German
umlauts, French accents, Turkish dotless i and Russian correctly with nothing
to download. The upstream face is 431 KB and mostly glyphs no locale in
SS_LANGUAGES uses, so it is subset here to 59 KB -- small enough to embed as a
byte array without making the build unpleasant.

Two things this script exists to keep honest:

  * The subset is a "Modified Version" under the SIL Open Font License, and
    Source Sans 3 reserves the font name "Source". Every name-table record
    carrying it therefore has to be rewritten, not just the filename. That is
    what most of the code below does.
  * The committed .ttf is a build artifact. This is its generator, in the same
    spirit as tools/codegen/ -- the output is committed so the build needs no
    network, and this script is how you check that it is what it claims to be.

The CJK faces are NOT built here: they are 16 MB each, they are fetched at
runtime (src/app/gui/FontCache.cpp), and they are shipped unmodified under
their own names, so nothing needs rewriting.

Usage:
    pip install fonttools
    python3 tools/make_ui_font.py            # downloads Source Sans 3, subsets
    python3 tools/make_ui_font.py --check    # verify the committed file matches
"""

import argparse
import hashlib
import io
import pathlib
import subprocess
import sys
import tempfile
import urllib.request
import zipfile

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

# The new identity. "Source" is a Reserved Font Name; a Modified Version may
# not carry it. See the OFL text shipped beside the font.
FAMILY = "Spirula UI"
STYLE = "Regular"
FULL_NAME = "Spirula UI"
PS_NAME = "SpirulaUI-Regular"
VERSION_STR = "Version 3.052;Spirula subset"
DESCRIPTION = (
    "Subset of Source Sans 3 (Adobe, SIL OFL 1.1) covering the Latin and "
    "Cyrillic ranges Spirula Studio's user interface uses. Renamed as the "
    "OFL requires of a Modified Version: 'Source' is a Reserved Font Name."
)

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "assets" / "fonts" / "SpirulaUI-Regular.ttf"


def fetch_source_sans() -> bytes:
    print(f"downloading {SOURCE_SANS_URL}")
    with urllib.request.urlopen(SOURCE_SANS_URL, timeout=120) as r:
        blob = r.read()
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        return z.read(SOURCE_SANS_MEMBER)


def rename(font) -> None:
    """Strip the Reserved Font Name from every name-table record."""
    name = font["name"]
    # nameID -> new value. 0 (copyright) and 13 (licence) must survive
    # unchanged; 7 is the trademark notice, which belongs to Adobe and must
    # not be reused, so it goes.
    replacements = {
        1: FAMILY,
        2: STYLE,
        3: f"{PS_NAME};{VERSION_STR}",
        4: FULL_NAME,
        5: VERSION_STR,
        6: PS_NAME,
        10: DESCRIPTION,
        16: FAMILY,
        17: STYLE,
    }
    for record in list(name.names):
        if record.nameID in (7, 8, 9, 11, 12, 19, 20, 21, 22, 25):
            name.names.remove(record)
            continue
        new = replacements.get(record.nameID)
        if new is not None:
            record.string = new
    # A leftover "Source" anywhere would defeat the whole exercise.
    for record in name.names:
        if record.nameID in (0, 13):
            continue          # copyright + licence: quoted verbatim, as required
        if "Source" in str(record):
            raise SystemExit(f"nameID {record.nameID} still says 'Source': {record}")


def build() -> bytes:
    from fontTools.ttLib import TTFont

    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        src = tmp / "in.ttf"
        dst = tmp / "out.ttf"
        src.write_bytes(fetch_source_sans())
        subprocess.run(
            [sys.executable, "-m", "fontTools.subset", str(src),
             f"--unicodes={UNICODES}",
             f"--layout-features={LAYOUT_FEATURES}",
             "--no-hinting", "--desubroutinize", "--drop-tables+=DSIG",
             f"--output-file={dst}"],
            check=True)
        # recalcTimestamp=False, or head.modified is stamped with the wall
        # clock on save and the committed file differs on every run -- which
        # would make --check useless.
        font = TTFont(dst, recalcTimestamp=False)
        rename(font)
        out = tmp / "final.ttf"
        font.save(out)
        return out.read_bytes()


def coverage_report(blob: bytes) -> None:
    from fontTools.ttLib import TTFont

    font = TTFont(io.BytesIO(blob))
    cmap = font.getBestCmap()
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
    print(f"  {'glyphs':12s} {len(font.getGlyphOrder())}, "
          f"{len(cmap)} mapped codepoints, {len(blob)} bytes")
    if bad:
        raise SystemExit("coverage check failed")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="rebuild and compare against the committed file")
    args = ap.parse_args()

    blob = build()
    coverage_report(blob)
    digest = hashlib.sha256(blob).hexdigest()

    if args.check:
        if not OUT.exists():
            raise SystemExit(f"{OUT} does not exist")
        have = hashlib.sha256(OUT.read_bytes()).hexdigest()
        if have != digest:
            raise SystemExit(f"{OUT} differs from a fresh build\n"
                             f"  committed {have}\n  rebuilt   {digest}")
        print(f"OK: {OUT.relative_to(ROOT)} matches ({digest[:16]}...)")
        return

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(blob)
    print(f"wrote {OUT.relative_to(ROOT)}  sha256 {digest}")


if __name__ == "__main__":
    main()
