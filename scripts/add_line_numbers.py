# -*- coding: utf-8 -*-
"""Add continuous line numbers to the left margin of the submission PDF.
Usage: python add_line_numbers.py <in.pdf> [out.pdf]
Writes in place (temp + replace) when out.pdf omitted.
"""
import sys, os, shutil

import pymupdf  # PyMuPDF

FS = 6.5            # line-number font size
GAP = 2.0           # gap between number and text block
MARGIN_X = 14.0     # x position of number (from left page edge)
COLOR = (0.35, 0.35, 0.35)
SKIP_TOP = 30.0     # ignore text rows starting above this y (running heads)
SKIP_BOTTOM = 40.0  # ignore rows starting below page_height - this (footers)


def page_lines(page):
    """Text lines in reading order with bboxes."""
    out = []
    d = page.get_text("dict")
    for block in d.get("blocks", []):
        if block.get("type") != 0:
            continue
        for line in block.get("lines", []):
            bbox = line.get("bbox")
            if not bbox:
                continue
            x0, y0, x1, y1 = bbox
            # skip running heads / footers / margin-overlapping rows
            if y0 < SKIP_TOP or y0 > page.rect.height - SKIP_BOTTOM:
                continue
            if x0 < MARGIN_X + FS + GAP:   # row starts too far left
                continue
            txt = "".join(s.get("text", "") for s in line.get("spans", []))
            if not txt.strip():
                continue
            out.append((y0, y1, x0))
    out.sort(key=lambda r: (round(r[0], 1), r[2]))
    # merge visually identical rows (same y) -> keep first
    merged, last_y = [], None
    for y0, y1, x0 in out:
        if last_y is not None and abs(y0 - last_y) < 1.5:
            continue
        merged.append((y0, y1, x0))
        last_y = y0
    return merged


def main():
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else src
    doc = pymupdf.open(src)
    n = 0
    for page in doc:
        for y0, y1, x0 in page_lines(page):
            n += 1
            ty = y0 + FS * 0.9
            page.insert_text((MARGIN_X, ty), str(n), fontsize=FS,
                             color=COLOR, fontname="helv")
    tmp = dst + ".tmp"
    doc.save(tmp, garbage=3, deflate=True)
    doc.close()
    if os.path.abspath(tmp) != os.path.abspath(dst):
        shutil.move(tmp, dst)
    else:
        os.replace(tmp, dst)
    print(f"OK: numbered {n} lines across document -> {dst}")


if __name__ == "__main__":
    main()
