# -*- coding: utf-8 -*-
"""Add native continuous line numbering to a .docx (w:lnNumType in sectPr).
Usage: python add_docx_line_numbers.py <in.docx>
In-place; safe to re-run (replaces existing lnNumType).
"""
import sys

from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement


def main():
    path = sys.argv[1]
    doc = Document(path)
    count = 0
    for section in doc.sections:
        sectPr = section._sectPr
        for old in sectPr.findall(qn("w:lnNumType")):
            sectPr.remove(old)
        ln = OxmlElement("w:lnNumType")
        ln.set(qn("w:countBy"), "1")
        ln.set(qn("w:restart"), "continuous")
        ln.set(qn("w:distance"), "360")
        sectPr.append(ln)
        count += 1
    doc.save(path)
    print(f"OK: continuous line numbering set on {count} section(s) -> {path}")


if __name__ == "__main__":
    main()
