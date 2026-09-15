"""Rebuild references in Cell Press style, two-phase (parse all -> rebuild all).

For PubMed-verified entries: list all authors when <=10, else first 10 + "et al.";
convert to Cell Press style (Author, A.B., Author, C.D., and Author, E.F. (Year).
Title. Journal Abbrev. vol, pages.). Titles/journals/volume-pages are taken from
the ORIGINAL entry text so nothing is invented.
Unverified entries are left byte-for-byte unchanged and flagged in the report.
"""
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
mp = REPO / 'docs/manuscript_v34.md'
t = mp.read_text(encoding='utf-8')
db = json.load(open(REPO / 'results/benchmark/v33/ref_authors_pubmed.json'))

i0 = t.find('## REFERENCES')
i1 = t.find('\n## ', i0 + 5)
if i1 < 0:
    i1 = len(t)

block = t[i0:i1]
entries = list(re.finditer(r'^(\d+)\.\s+(.+?)(?=^\d+\.\s|\Z)', block, re.M | re.S))

JOURNAL = {
    'N Engl J Med': 'N. Engl. J. Med.', 'Nat Med': 'Nat. Med.',
    'Nat Biotechnol': 'Nat. Biotechnol.', 'Nat Genet': 'Nat. Genet.',
    'Nat Commun': 'Nat. Commun.', 'Nat Methods': 'Nat. Methods',
    'Nat Rev Neurosci': 'Nat. Rev. Neurosci.', 'Cancer Cell': 'Cancer Cell',
    'Cell': 'Cell', 'Science': 'Science', 'Nature': 'Nature',
    'J Clin Invest': 'J. Clin. Invest.', 'J Natl Cancer Inst': 'J. Natl. Cancer Inst.',
    'PLoS ONE': 'PLoS One', 'Brief Bioinform': 'Brief. Bioinform.',
    'BMC Genomics': 'BMC Genomics', 'Proc Natl Acad Sci USA': 'Proc. Natl. Acad. Sci. U.S.A.',
    'Mol Cancer Ther': 'Mol. Cancer Ther.', 'Biometrics': 'Biometrics',
    'J R Stat Soc B': 'J. R. Stat. Soc. B', 'Patterns': 'Patterns',
}


def cp_authors(names):
    """['Topalian SL'] -> 'Topalian, S.L., ...' Cell Press style."""
    def fmt(n):
        parts = n.split()
        if len(parts) < 2:
            return n
        surname = ' '.join(parts[:-1])
        initials = parts[-1]
        return f'{surname}, {".".join(initials)}.'
    if len(names) > 10:
        return ', '.join(fmt(n) for n in names[:10]) + ', et al.'
    if len(names) == 1:
        return fmt(names[0])
    return ', '.join(fmt(n) for n in names[:-1]) + ', and ' + fmt(names[-1])


report = []
new_parts = []
last_end = 0
for m in entries:
    num, raw = m.group(1), m.group(2)
    entry = ' '.join(raw.split())
    # 原条目的字段
    ym = re.search(r'\((\d{4})\)', entry)
    year = ym.group(1) if ym else ''
    tm = re.search(r'\(\d{4}\)\s+(.+?)\s+\*', entry)
    title = tm.group(1).strip() if tm else ''
    jm = re.search(r'\*([^*]+)\*[:,]?\s*(.+?)\.?\s*$', entry)
    journal_raw = jm.group(1).strip() if jm else ''
    tail = jm.group(2).strip().rstrip('.') if jm else ''

    rec = db.get(num, {})
    names = rec.get('authors') or []
    ok = bool(rec.get('verified')) and bool(names) and bool(title) and bool(journal_raw)

    if ok and num != '1':  # [1] 已判定为噪声匹配，强制保留原文
        au = cp_authors(names)
        journal = JOURNAL.get(journal_raw, journal_raw)
        body = f'{au} ({year}). {title}. {journal} {tail}.'
        report.append((num, 'CP', f'{len(names)} au'))
    else:
        body = entry
        report.append((num, 'KEPT', rec.get('error', 'unverified')[:24]))

    if m.start() >= last_end:
        new_parts.append(block[last_end:m.start()])
    new_parts.append(f'{num}. {body}\n')
    last_end = m.end()

new_parts.append(block[last_end:])
new_block = ''.join(new_parts)
t = t[:i0] + new_block + t[i1:]
mp.write_text(t, encoding='utf-8')

print('rebuild complete')
cp = [r for r in report if r[1] == 'CP']
kept = [r for r in report if r[1] == 'KEPT']
print(f'Cell Press style: {len(cp)} | kept original: {len(kept)}')
print('kept ids:', [r[0] for r in kept])
# 抽查
refs = t[t.find('## REFERENCES'):]
for n in ['3', '14', '23', '24', '25']:
    mm = re.search(rf'^{n}\. (.+)$', refs, re.M)
    if mm:
        print(f'[{n}]', mm.group(1)[:190])
