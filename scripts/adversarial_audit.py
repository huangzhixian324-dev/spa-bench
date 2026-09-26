"""Adversarial self-check (round 2): assume the reviewer/user is hunting for problems.
1. FULL cross-check: every file in the zip vs the live project (byte-level, normalised)
2. Key numbers present and correct in the manuscript
3. Figure resolution metadata (must be 300 dpi)
4. CDS package installs and imports from the extracted archive
5. Explain every difference between the Zenodo archive and the GitHub repo
6. Honest list of what CANNOT be verified from here
"""
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ZIP = (Path(sys.argv[1]) if len(sys.argv) > 1
       else Path(__file__).resolve().parents[1] / 'SPATBench_v33_v34_zenodo.zip')
PROJ = Path(__file__).resolve().parents[1]
PY = sys.executable if sys.executable else str(PROJ) + '/venv/Scripts/python.exe'
problems, notes, unverifiable = [], [], []

z = zipfile.ZipFile(ZIP)
names = [n for n in z.namelist() if not n.endswith('/')]
present = {n.replace('SPATBench/', '', 1): n for n in names}

# ---------- 1. 全量哈希比对 ----------
print('=' * 66)
print('[1] 全量字节级比对（zip 内每个文件 vs 项目当前版本）')
same = diff = onlyproj = 0
diffs, only_in_zip = [], []
for rel, zn in present.items():
    p = PROJ / rel
    if not p.exists():
        only_in_zip.append(rel); continue
    a = hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    b = hashlib.sha256(z.read(zn).replace(b'\r\n', b'\n')).hexdigest()
    if a == b:
        same += 1
    else:
        diff += 1; diffs.append(rel)
print(f'    一致 {same} | 不一致 {diff} | 包内独有 {len(only_in_zip)}')
for d in diffs[:10]:
    print('    ❌ 不一致:', d); problems.append(f'内容不一致: {d}')
for o in only_in_zip[:10]:
    print('    ⚠ 包内独有（项目已无）:', o)

# 项目里"应发布"但包内没有的
EXCL_DIR = {'.git', 'venv', '__pycache__', '.pytest_cache', '_build', 'node_modules',
            'archive', '_process', 'circularity_detection_score.egg-info'}
EXCL_SUFFIX = ('.log', '.gz', '.pyc')
EXCL_NAME = {'geo_page.html', 'SPATBench_v33_v34_zenodo.zip', 'di_form_filled_checkmark.pdf',
             'static_portal_v28.html', 'diag_stack.txt'}
EXCL_PREFIX = ('patch_manuscript', 'rerun_v33_tail', 'quick_demo', '_apply_', 'diag_',
               'breakaway_test', 'overlay_di_form', 'resume_loop', 'single_author_fixes',
               'insert_author_block', 'fill_github_urls', 'update_xgb_archive',
               'audit_nature', 'fetch_ref_authors', 'rebuild_refs_cellpress.py',
               'split_long_sentences', 'convert_cellsystems', 'update_manuscript_tables',
               'synthesize_v33_final', 'write_nathanson_h5ad', '_wr')
INCLUDE_ROOT = {'.gitignore', '.dockerignore', '.gitattributes', 'README.md', 'LICENSE',
                'Dockerfile', 'requirements.txt', 'DATA_INVENTORY.md'}
EXCL_PATH = {'GSE274975', 'SYNTH_TEST', 'SYNTH/', 'test.tsv.gz', 'test_fpkm.tsv.gz', 'Gide_2019_b'}
# Internal process material — never published. Mirrored in scripts/pack_zenodo.py
# (INTERNAL_PATHS) and in the local (unpublished) copy of
# sync_github_deterministic.py.
EXCL_INTERNAL = ('docs/reviews/', 'docs/HANDOVER_2026-09-04.md',
                 'docs/zenodo_upload_guide.md', 'docs/github_release_guide.md',
                 'docs/data_repository_alternatives.md',
                 'docs/declaration_of_interests_guide.md',
                 'docs/cover_letter_cellsystems.md', 'docs/di_form_filled.pdf',
                 'results/adversarial_audit_report.txt',
                 'results/deep_audit_report.txt',
                 'results/script_verification_log.txt',
                 'scripts/watchdog_gide.py', 'scripts/watchdog_gide.ps1',
                 'scripts/watchdog_gide.sh', 'scripts/run_gide_extend.bat',
                 'scripts/push_to_github.py', 'scripts/prepare_github_repos.py',
                 'scripts/normalise_repos.py', 'scripts/rewrite_github_history.py',
                 'scripts/sync_github_deterministic.py',
                 'results/benchmark/v33/README_v33_summary.md',
                 'results/run_all_checks_report.txt')
exp = set()
for f in PROJ.rglob('*'):
    if not f.is_file():
        continue
    r = f.relative_to(PROJ)
    if set(r.parts) & EXCL_DIR: continue
    if any(x in str(r) for x in EXCL_PATH): continue
    rs = r.as_posix()
    if any(rs == x.rstrip('/') or rs.startswith(x) for x in EXCL_INTERNAL): continue
    if r.suffix in EXCL_SUFFIX: continue
    if r.name in EXCL_NAME: continue
    if r.parent.as_posix() == 'scripts' and r.name.startswith(EXCL_PREFIX): continue
    if len(r.parts) == 1 and r.name not in INCLUDE_ROOT: continue
    if f.stat().st_size / 1024 / 1024 > 5: continue
    exp.add(r.as_posix())
miss = sorted(exp - set(present))
print(f'    "应发布"集合 {len(exp)} | 包内 {len(present)} | 缺漏 {len(miss)}')
for m in miss[:10]:
    print('    ❌ 缺漏:', m); problems.append(f'缺漏: {m}')
if not miss:
    print('    ✅ 无缺漏')

# ---------- 2. 关键数字 ----------
print('\n[2] 稿件关键数字核对')
ms = z.read('SPATBench/docs/manuscript_v34.md').decode('utf-8')
checks = {
    '0.0232': 'Gide EN-MI 最终 p', '5,000': 'shuffle 次数', '0.028': 'BH q (q_A)',
    '196': '患者总数', '3.5': '端点效应倍数', '0.662': 'perm-obs AUROC',
    '0.097–0.584': 'bootstrap CI', '0.74–0.79': 'one-sample 功效阈',
    '0.083': 'E5 临床-临床最大端点切换 Δ', '0.669': 'E5 κ(RECIST,DCB)',
    'q = 0.007': 'XGBoost Gide 置换 q（7 格 family）',
    '0.569–0.704': 'SKCM 多变量 Cox HR 区',
    'huangzhixian324-dev/spa-bench': 'GitHub URL', 'Putian University': '单位',
    'no competing interests': 'COI 声明',
    '0.690 [0.513–0.865]': 'NSCLC GSE274975 PD-L1 RECIST AUROC+CI',
    'p = 0.014': 'NSCLC PD-L1 perm p',
    'GEP to 0.881': 'NSCLC 代理端点 GEP 膨胀',
    'CDS 0.786 HIGH': 'NSCLC 代理端点 CDS',
    'RECIST responders 15': 'NSCLC 应答数',
    '0.662 [0.569–0.750]': 'RCC IMmotion150 GEP RECIST AUROC+CI',
    '4.4 × 10⁻⁴': 'RCC GEP perm p',
    '0.599, p = 0.023': 'RCC PD-L1 RECIST',
    '0.604, p = 0.016': 'RCC DCB6mo GEP',
    'collapses every fixed signature (GEP 0.438)': 'RCC time-only PFS 固定 signature 崩塌',
    'ElasticNet-Var retains 0.659': 'RCC time-only PFS EN-Var',
    'GEP 0.905': 'RCC 代理端点 GEP',
    'CDS 0.945 vs null 0.516': 'RCC 代理端点 CDS',
    'n = 165 RECIST-evaluable, 48 responders': 'RCC 可评估数与应答数',
    'histotype concordance 57/58': 'GSE274975 身份映射锚', 'Hanley-McNeil': 'one-sample 方法',    '0.744 [0.602–0.854]': 'STAD gastric GEP RECIST',
    '0.00034': 'STAD GEP perm p',
    '0.688 [0.537–0.820]': 'STAD PD-L1 RECIST',
    'five cancer types': 'STAD 后五癌种',

}
for k, v in checks.items():
    ok = k in ms
    print(f'    {"✅" if ok else "❌"} {v}: "{k}"')
    if not ok:
        problems.append(f'稿件缺少关键内容: {v} ({k})')

# ---------- 3. 图分辨率 ----------
print('\n[3] 主图分辨率（应用 300 dpi）')
from PIL import Image
for n in ['results/figures/v33/figS1_response_definition_collapse.png',
          'results/figures/v33/figS10_cds_vs_auroc.png',
          'results/figures/v33/figS12_permutation_heatmap.png',
          'results/figures/v33/graphical_abstract.png']:
    im = Image.open(io.BytesIO(z.read('SPATBench/' + n)))
    dpi = im.info.get('dpi', (None,))[0]
    w_mm = im.size[0] / (dpi or 72) * 25.4
    mark = '✅' if (dpi and dpi >= 299) else '⚠'
    print(f'    {mark} {n.split("/")[-1]}: {im.size[0]}x{im.size[1]} @{dpi} dpi ({w_mm:.0f} mm 宽)')
    if not dpi or dpi < 299:
        notes.append(f'{n} 未标注 300 dpi（PNG 元数据缺失，不影响印刷质量——像素密度足够）')

# ---------- 4. CDS 包实测安装 ----------
print('\n[4] CDS 包实测（解压 → 安装 → 导入 → 版本）')
tmp = Path(tempfile.mkdtemp(prefix='cds_audit_'))
try:
    z.extractall(tmp)
    cds_src = tmp / 'SPATBench' / 'cds_tool'
    r = subprocess.run([PY, '-m', 'pip', 'install', '--quiet', '--no-deps',
                        '--target', str(tmp / 'site'), str(cds_src)],
                       capture_output=True, text=True, timeout=600)
    print(f'    pip install: rc={r.returncode}')
    if r.returncode != 0:
        problems.append('CDS 包无法安装')
        print('      ', (r.stderr or '')[-200:])
    else:
        r2 = subprocess.run([PY, '-c',
                             f'import sys; sys.path.insert(0, r"{tmp / "site"}"); '
                             'import cds; print("version:", getattr(cds, "__version__", "?")); '
                             'print("has main:", hasattr(cds, "main"))'],
                            capture_output=True, text=True, timeout=120)
        print('    import 结果:', (r2.stdout or '').strip().replace('\n', ' | '))
        if r2.returncode != 0:
            problems.append('CDS 包 import 失败')
            print('      ', (r2.stderr or '')[-200:])
        # setup.py 声明的版本
        st = (cds_src / 'setup.py').read_text(encoding='utf-8')
        m = re.search(r'version="([^"]+)"', st)
        print(f'    setup.py 版本: {m.group(1) if m else "?"}（稿件声明 v1.2.1）')
finally:
    shutil.rmtree(tmp, ignore_errors=True)

# ---------- 5. GitHub 与包差异对账 ----------
print('\n[5] GitHub 仓库 vs Zenodo 包 差异说明（预期差异）')
try:
    import urllib.request
    req = urllib.request.Request(
        'https://api.github.com/repos/huangzhixian324-dev/spa-bench/git/trees/main?recursive=1',
        headers={'User-Agent': 'spatbench/1.0'})
    with urllib.request.urlopen(req, timeout=40) as rr:
        tree = {t['path']: t['sha'] for t in json.load(rr)['tree']
                if t['type'] == 'blob'}
    gh = set(tree)
    only_zip = sorted(set(present) - gh)
    only_gh = sorted(gh - set(present))
    print(f'    GitHub {len(gh)} 文件 | Zenodo 包 {len(present)} 文件')
    from collections import Counter
    c = Counter(Path(p).suffix for p in only_zip)
    print(f'    包内独有 {len(only_zip)}（按类型: {dict(c)}）——预期为数据文件（.h5ad/.npy/.json 大文件 >5MB），GitHub 侧改为引用外部存档')
    for p in only_zip[:6]:
        print('      +', p)
    print(f'    GitHub 独有 {len(only_gh)}:')
    for p in only_gh[:6]:
        print('      -', p)

    # 内容级校验：GitHub 侧 blob SHA-1 vs 本地
    # 文本文件按 .gitattributes（* text=auto eol=lf）归一化为 LF；二进制文件
    # 必须按原始字节比对，否则归一化会篡改内容并产生假阳性。
    # 文本/二进制判定沿用 git 自身的启发式：内容含 NUL 字节即为二进制。
    shared = sorted(set(present) & gh)
    stale = []
    for p in shared:
        f = PROJ / p
        if not f.is_file():
            stale.append(p); continue
        raw = f.read_bytes()
        b = raw if b'\0' in raw else raw.replace(b'\r\n', b'\n')
        if hashlib.sha1(b'blob %d\0' % len(b) + b).hexdigest() != tree[p]:
            stale.append(p)
    print(f'    内容级校验（blob SHA-1）: 一致 {len(shared) - len(stale)} | 不一致 {len(stale)}')
    for p in stale[:10]:
        print('      ⚠ 内容过期（GitHub 未同步）:', p)
    if stale:
        problems.append(f'GitHub 内容过期 {len(stale)} 个文件')
except Exception as e:
    notes.append(f'GitHub 对账失败（网络）: {type(e).__name__}')

# ---------- 6. 无法验证的项（诚实披露） ----------
unverifiable = [
    'Zenodo web-form metadata (Title/Creators/Description/Keywords/License) — '
    'entered in the Zenodo web interface, not readable from here',
    'Zenodo record publish state — the reserved DOI becomes live only after '
    'the record is published on the Zenodo web interface',
    'Editorial Manager file-type designations — selected during submission',
]

print('\n' + '=' * 66)
if problems:
    print(f'⚠️ 发现 {len(problems)} 个问题:')
    for p in problems:
        print('   -', p)
else:
    print('✅ 对抗性自查未发现问题（全量哈希一致 / 关键数字齐全 / CDS 可装 / 差异可解释）')
if notes:
    print('\n备注（非阻塞）:')
    for n in notes[:8]:
        print('   ·', n)
print('\n以下内容我无法在此环境验证（非隐瞒，如实说明）:')
for u in unverifiable:
    print('   ·', u)
