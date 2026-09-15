"""Apply the Zenodo DOI to the manuscript (one command, ready for when the DOI arrives).
Usage: python apply_doi.py 10.5281/zenodo.XXXXXXX
Does: 1) fill Data-and-code availability, 2) add a Key Resources Table row,
3) enhance the Data bullet with the archive DOI, 4) regenerate PDF + Word,
5) verify every occurrence.
"""
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MP = REPO / 'docs' / 'manuscript_v34_nar.md'
PY = sys.executable if sys.executable else str(R) + '/venv/Scripts/python.exe'


def main(doi: str):
    doi = doi.strip()
    if not re.fullmatch(r'10\.\d{4,9}/\S+', doi):
        print('DOI 格式可疑:', doi, '（应形如 10.5281/zenodo.1234567）')
    url = 'https://doi.org/' + doi if not doi.startswith('http') else doi
    t = MP.read_text(encoding='utf-8')
    n = 0

    # 1) Data and code availability — fill the placeholder
    old = 'DOI: [ZENODO DOI — to be inserted at acceptance]'
    new = f'DOI: {doi}'
    if old in t:
        t = t.replace(old, new); n += 1; print('[1] Data availability DOI 已填')
    else:
        print('[1] 占位未找到（可能已填）')

    # 2) Data bullet — mention the archived preprocessed inputs
    old2 = ('- Data: All raw cohort data are publicly available from GEO '
            '(GSE78220, GSE91061, GSE100797, GSE135222) and cBioPortal '
            '(mel_iatlas_gide_2019, blca_iatlas_imvigor210_2017, mel_iatlas_liu_2019).')
    new2 = (old2 + ' Preprocessed expression matrices, all per-cell permutation '
            'outputs and figures are archived in the versioned Zenodo record '
            f'(DOI: {doi}).')
    if old2 in t:
        t = t.replace(old2, new2); n += 1; print('[2] Data 行已补存档说明')

    # 3) Key Resources Table — add a Zenodo row
    krt_anchor = ('| circularity-detection-score | This paper; GitHub '
                  '(https://github.com/huangzhixian324-dev/cds) | v1.2.1 |')
    if krt_anchor in t and 'SPATBench archive (Zenodo)' not in t:
        t = t.replace(krt_anchor, krt_anchor +
                      f'\n| SPATBench archive (Zenodo) | This paper | {doi} |')
        n += 1
        print('[3] KRT 已加 Zenodo 行')

    MP.write_text(t, encoding='utf-8')
    print(f'共 {n} 处修改')

    # 4) regenerate
    for s in ['build_submission.py', 'md_to_docx_cellpress.py']:
        r = subprocess.run([PY, str(REPO / 'scripts' / s)], cwd=str(REPO),
                           capture_output=True, text=True, timeout=900)
        print(f'   {s}: rc={r.returncode}')

    # 5) verify
    from pypdf import PdfReader
    pdf = REPO / 'docs' / 'SPATBench_CellSystems_submission_v35.pdf'
    txt = '\n'.join((p.extract_text() or '') for p in PdfReader(pdf).pages)
    print(f'\nPDF 中 DOI 出现次数: {txt.count(doi)}  （页面 {len(PdfReader(pdf).pages)}）')
    t2 = MP.read_text(encoding='utf-8')
    print(f'Markdown 中 DOI 出现次数: {t2.count(doi)}')
    print(f'\n✅ 完成。DOI = {doi}')
    print(f'   引用链接: {url}')
    (REPO / 'results' / 'zenodo_doi.txt').write_text(doi, encoding='utf-8')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('用法: python apply_doi.py 10.5281/zenodo.XXXXXXX')
    else:
        main(sys.argv[1])
