"""Assemble the NMI Analysis submission PDFs from the NMI manuscript markdown.

Parallel to scripts/build_submission.py (Cell Systems) but:
  - sources: docs/manuscript_nmi_v1.md + docs/supplementary_nmi_v1.md
  - outputs: docs/SPATBench_NMI_analysis_v1.pdf / _supplementary_v1.pdf
  - no graphical abstract insertion (not required by NMI)
Outputs the same reviewer-readable typeset PDFs with figures above legends.
"""
import re
import subprocess
from pathlib import Path
import markdown

REPO = Path(__file__).resolve().parents[1]
MD = REPO / 'docs' / 'manuscript_nmi_v1.md'
FIGD = REPO / 'results' / 'figures' / 'v33'
HTML_OUT = REPO / 'docs' / '_build' / 'submission_nmi.html'
HTML_OUT.parent.mkdir(parents=True, exist_ok=True)
PDF_OUT = REPO / 'docs' / 'SPATBench_NMI_analysis_v1.pdf'
EDGE = r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'

t = MD.read_text(encoding='utf-8')

# figures above their legends (legend anchors identical to v34 back matter)
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

body = markdown.markdown(t, extensions=['tables', 'fenced_code'])
body = re.sub(r'\[(\d+(?:[,&\u2013\-]\s?\d+)*)\]', r'<sup>[\1]</sup>', body)

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

r = subprocess.run(
    [EDGE, '--headless', '--disable-gpu',
     '--print-to-pdf=' + str(PDF_OUT),
     '--no-pdf-header-footer',
     'file:///' + str(HTML_OUT).replace('\\', '/')],
    capture_output=True, text=True, timeout=120)
print('edge rc=', r.returncode)
print('PDF written:', PDF_OUT, PDF_OUT.stat().st_size // 1024, 'KB')

# supplementary (NMI version, includes S31)
SUP_MD = REPO / 'docs' / 'supplementary_nmi_v1.md'
SUP_HTML = REPO / 'docs' / '_build' / 'supplementary_nmi.html'
SUP_PDF = REPO / 'docs' / 'SPATBench_NMI_supplementary_v1.pdf'

sup_lines = SUP_MD.read_text(encoding='utf-8').split('\n')


def _fig_for(head_idx):
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
