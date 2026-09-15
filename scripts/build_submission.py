"""Assemble the Cell Systems submission PDFs from manuscript markdown.
- Main text: inserts graphical abstract + Figures 1/2/3 above their legends
- Supplementary: inserts Figures S1-S12 above their captions
- Renders Markdown (tables extension) to styled single-column HTML
- Edge headless converts HTML -> PDF
Outputs:
  docs/SPATBench_CellSystems_submission_v35.pdf      (main text)
  docs/SPATBench_CellSystems_supplementary_v35.pdf   (supplementary material)
"""
import re
import subprocess
from pathlib import Path
import markdown

REPO = Path(__file__).resolve().parents[1]
MD = REPO / 'docs' / 'manuscript_v34.md'
FIGD = REPO / 'results' / 'figures' / 'v33'
HTML_OUT = REPO / 'docs' / '_build' / 'submission.html'
HTML_OUT.parent.mkdir(parents=True, exist_ok=True)
PDF_OUT = REPO / 'docs' / 'SPATBench_CellSystems_submission_v35.pdf'
EDGE = r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'

t = (REPO / 'docs' / 'manuscript_v34.md').read_text(encoding='utf-8')

# 1. graphical abstract after the H1 title line
ga = FIGD / 'graphical_abstract.png'
# 图形摘要插在作者块（含 affiliations 行）之后
anchor = '**Author affiliations:** ^1 [same as above]'
if anchor in t:
    t = t.replace(anchor, anchor + '\n\n![](file:///'
                  + str(ga).replace('\\', '/') + ')\n\n'
                  '**Graphical abstract.** CDS workflow (left) and the three main '
                  'findings (right); bottom: high-resolution validation on Gide 2019.',
                  1)
else:
    t = re.sub(r'(# SPATBench:.*\n)',
               r'\1\n\n![](file:///' + str(ga).replace('\\', '/') + ')\n\n'
               '**Graphical abstract.** CDS workflow (left) and the three main '
               'findings (right); bottom: high-resolution validation on Gide 2019.\n',
               t, count=1)

# 2. figures above their legends
figs = [
    ('**Figure 1. Response-definition collapse',
     FIGD / 'figS1_response_definition_collapse.png'),
    ('**Figure 2. CDS versus observed trainable-ML performance.',
     FIGD / 'figS10_cds_vs_auroc.png'),
    ('**Figure 3. Complete permutation test matrix.',
     FIGD / 'figS12_permutation_heatmap.png'),
]
for anchor, img in figs:
    assert anchor in t, anchor
    t = t.replace(anchor,
                  '![](file:///' + str(img).replace('\\', '/') + ')\n\n' + anchor,
                  1)

# 3. markdown -> html
body = markdown.markdown(t, extensions=['tables', 'fenced_code'])

# Cell Systems citation style: in-text [n] as superscript (outside refs list)
import re as _re
body = _re.sub(r'\[(\d+(?:[,&]\s?\d+)*)\]', r'<sup>[\1]</sup>', body)
# keep the References section numbering as-is (it uses "1." list markers)

CSS = """
@page { size: A4; margin: 22mm 20mm; }
body { font-family: 'Times New Roman', 'STIX Two Text', serif;
       font-size: 11pt; line-height: 1.45; color: #111; }
h1 { font-size: 17pt; text-align: left; margin: 0 0 6mm; }
h2 { font-size: 13.5pt; border-bottom: 1px solid #999; padding-bottom: 2px;
     margin-top: 9mm; }
p  { margin: 2mm 0; text-align: justify; }
img { max-width: 100%; display: block; margin: 4mm auto 2mm; }
table { border-collapse: collapse; width: 100%; font-size: 8.5pt;
        margin: 3mm 0; }
th, td { border: 1px solid #666; padding: 2.5px 5px; text-align: left; }
th { background: #eee; }
code { font-family: Consolas, monospace; font-size: 9pt; }
strong { font-weight: bold; }
em { font-style: italic; }
hr { border: none; border-top: 1px solid #bbb; margin: 6mm 0; }
ol, ul { margin: 2mm 0; padding-left: 8mm; }
"""

html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<style>{CSS}</style></head><body>{body}</body></html>"""
HTML_OUT.write_text(html, encoding='utf-8')
print('HTML written:', HTML_OUT)

# 4. Edge headless -> PDF
r = subprocess.run(
    [EDGE, '--headless', '--disable-gpu',
     '--print-to-pdf=' + str(PDF_OUT),
     '--no-pdf-header-footer',
     'file:///' + str(HTML_OUT).replace('\\', '/')],
    capture_output=True, text=True, timeout=120)
print('edge rc=', r.returncode)
print('PDF written:', PDF_OUT, PDF_OUT.stat().st_size // 1024, 'KB')

# 5. supplementary material -> reviewer-readable PDF
#    (docs/supplementary_material_v33.md is the generator's source; reviewers
#     must receive a typeset document, not raw Markdown)
SUP_MD = REPO / 'docs' / 'supplementary_material_v33.md'
SUP_HTML = REPO / 'docs' / '_build' / 'supplementary.html'
SUP_PDF = REPO / 'docs' / 'SPATBench_CellSystems_supplementary_v35.pdf'

sup_lines = SUP_MD.read_text(encoding='utf-8').split('\n')


def _fig_for(head_idx):
    """Filename referenced by the caption block that follows a figure heading."""
    for ln in sup_lines[head_idx + 1:]:
        if ln.startswith('## '):
            break
        m = re.search(r'`(figS\d+[^`]*\.png)`', ln)
        if m:
            return m.group(1)
    return None


out_lines, inserted = [], 0
for i, ln in enumerate(sup_lines):
    out_lines.append(ln)
    if re.match(r'^### Figure S\d+\.', ln):
        fn = _fig_for(i)
        if fn and (FIGD / fn).exists():
            out_lines += ['', '![](file:///' + str(FIGD / fn).replace('\\', '/') + ')']
            inserted += 1
print('supplementary figures embedded:', inserted)
assert inserted >= 10, f'expected >=10 supplementary figures, got {inserted}'

sup_body = markdown.markdown('\n'.join(out_lines),
                             extensions=['tables', 'fenced_code'])
SUP_HTML.write_text(
    f'<!DOCTYPE html><html><head><meta charset="utf-8">'
    f'<style>{CSS}</style></head><body>{sup_body}</body></html>',
    encoding='utf-8')
print('HTML written:', SUP_HTML)

r = subprocess.run(
    [EDGE, '--headless', '--disable-gpu',
     '--print-to-pdf=' + str(SUP_PDF),
     '--no-pdf-header-footer',
     'file:///' + str(SUP_HTML).replace('\\', '/')],
    capture_output=True, text=True, timeout=180)
print('edge rc=', r.returncode)
print('PDF written:', SUP_PDF, SUP_PDF.stat().st_size // 1024, 'KB')
