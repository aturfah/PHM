#!/usr/bin/env python3
"""
build.py -- assemble the standalone PHM Explorer page.

template.html holds everything (styles + markup + app logic) but is a bare
fragment, not a full document, and its <script id="app-data"> placeholder
(__DATA_JSON__) still needs a real export to be useful out of the box.

This script does two things:
  1. Wraps the fragment in a proper <!doctype html><html><head>...<body>
     document (everything up to the first <div class="page"> becomes the
     <head>; everything from there on becomes the <body>).
  2. Substitutes __DATA_JSON__ with the contents of an export JSON file, so
     the page opens with real data already loaded instead of an empty shell.

Usage:
    python3 build.py [--data ../examples/two_moons_export.json] [--out ../index.html]

The generated file is completely standalone -- open it directly in a
browser, no server or build step required. Loading a different export at
runtime (via the "load a PHM export..." picker in the page itself) never
requires rebuilding; this script only controls which dataset ships as the
*default* when the page is opened cold.
"""
import argparse
import pathlib

HERE = pathlib.Path(__file__).resolve().parent

def build(template_path, data_path, out_path):
    template = template_path.read_text(encoding="utf-8")
    data = data_path.read_text(encoding="utf-8")

    if "__DATA_JSON__" not in template:
        raise SystemExit(f"{template_path} has no __DATA_JSON__ placeholder")

    marker = '<div class="page">'
    split_at = template.find(marker)
    if split_at == -1:
        raise SystemExit(f"{template_path} has no {marker!r} to split head/body on")

    head_fragment = template[:split_at]
    body_fragment = template[split_at:]

    # template.html already opens with its own <meta charset="utf-8"> (needed
    # when it's inlined as an Artifact fragment); drop it here since the
    # wrapper above already declares one, to avoid a duplicate.
    head_fragment = head_fragment.replace('<meta charset="utf-8">\n', "", 1)

    doc = (
        "<!doctype html>\n"
        '<html lang="en">\n'
        "<head>\n"
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"{head_fragment}"
        "</head>\n"
        "<body>\n"
        f"{body_fragment}\n"
        "</body>\n"
        "</html>\n"
    )
    doc = doc.replace("__DATA_JSON__", data)

    out_path.write_text(doc, encoding="utf-8")
    print(f"wrote {out_path} ({out_path.stat().st_size:,} bytes)")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--template", default=str(HERE / "template.html"))
    p.add_argument("--data", default=str(HERE.parent / "examples" / "two_moons_export.json"),
                   help="export JSON to embed as the page's default dataset")
    p.add_argument("--out", default=str(HERE.parent / "index.html"))
    args = p.parse_args()

    build(pathlib.Path(args.template), pathlib.Path(args.data), pathlib.Path(args.out))
