"""Deep audit of the archive (optional, slower).

Usage:
    python scripts/deep_audit_archive.py [path/to/SPATBench_v33_v34_zenodo.zip]

1. every JSON / PNG / H5AD member actually parsed (not sampled)
2. required-file checklist (missing = hard fail)
3. extract to a temp directory and RUN the pipeline self-test (validate.py)
4. spot-check byte-level consistency against the live repository
5. path / encoding sanity
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

REPO = Path(__file__).resolve().parents[1]
ZIP = (Path(sys.argv[1]) if len(sys.argv) > 1
       else REPO / 'SPATBench_v33_v34_zenodo.zip')

if not ZIP.exists():
    print(f'找不到归档文件: {ZIP}')
    print('用法: python scripts/deep_audit_archive.py '
          '[path/to/SPATBench_v33_v34_zenodo.zip]')
    print('（无参数时默认在仓库根目录查找）')
    sys.exit(2)

print(f'检查目标: {ZIP.name}  {ZIP.stat().st_size / 1024 / 1024:.1f} MB\n')
problems = []
notes = []

z = zipfile.ZipFile(ZIP)
names = [n for n in z.namelist() if not n.endswith('/')]
present = {n.replace('SPATBench/', '', 1) for n in names}
print(f'文件数: {len(names)}')


def inner(n):
    return n.replace('SPATBench/', '', 1)


# ---------- 1. 全量格式校验 ----------
print('\n[1] 全量格式校验')
cnt = {'json': 0, 'png': 0, 'h5ad': 0, 'npy': 0}
for n in names:
    low = n.lower()
    try:
        if low.endswith('.json'):
            json.loads(z.read(n)); cnt['json'] += 1
        elif low.endswith('.png'):
            from PIL import Image
            im = Image.open(io.BytesIO(z.read(n))); im.load(); cnt['png'] += 1
        elif low.endswith('.h5ad'):
            import h5py
            with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as t:
                t.write(z.read(n)); tp = t.name
            with h5py.File(tp, 'r') as f:
                assert 'X' in f or 'obs' in f
            os.unlink(tp); cnt['h5ad'] += 1
        elif low.endswith('.npy'):
            import numpy as np
            np.load(io.BytesIO(z.read(n)), allow_pickle=True); cnt['npy'] += 1
    except Exception as e:
        problems.append(f'格式损坏 {inner(n)}: {str(e)[:60]}')
for k, v in cnt.items():
    print(f'    {k}: {v} 个正常' + ('（有损坏，见结论）' if problems else ' ✅'))

# ---------- 2. 必备文件 ----------
print('\n[2] 必备文件清单')
required = [
    'docs/manuscript_v34_nar.md', 'docs/supplementary_material_v33.md',
    'docs/SPATBench_CellSystems_submission_v35.pdf',
    'docs/SPATBench_CellSystems_manuscript_v35.docx',
    'results/figures/v33/figS1_response_definition_collapse.png',
    'results/figures/v33/figS12_permutation_heatmap.png',
    'results/figures/v33/graphical_abstract.png',
    'results/benchmark/v33/benchmark_v33.json',
    'results/benchmark/v33/permutation_v33.json',
    'results/benchmark/v33/nc_cells/Gide_2019_cBio__ElasticNet_5000.json',
    'results/benchmark/v33/effect_size_bootstrap_v33.json',
    'results/benchmark/v33/power_onesample_v33.json',
    'results/benchmark/v33/xgb_exclusion_evidence_v33.json',
    'scripts/rerun_v33.py', 'scripts/validate.py', 'scripts/make_figures_v33.py',
    'scripts/build_submission.py', 'scripts/run_all_checks.py',
    'cds_tool/cds.py', 'cds_tool/setup.py', 'cds_tool/LICENSE',
    'README.md', 'LICENSE', 'Dockerfile', 'requirements.txt',
]
missing = [r for r in required if r not in present]
print(f'    清单 {len(required)} 项，缺失 {len(missing)}')
for m in missing:
    problems.append(f'缺失: {m}')
    print('    ❌', m)
if not missing:
    print('    ✅ 全部在位')

# ---------- 3. 解包实测 ----------
print('\n[3] 解包实测（运行包内 validate.py，慢）')
tmp = Path(tempfile.mkdtemp(prefix='spatbench_deep_'))
try:
    z.extractall(tmp)
    root = tmp / 'SPATBench'
    r = subprocess.run([sys.executable, str(root / 'scripts' / 'validate.py')],
                       cwd=str(root), capture_output=True, text=True, timeout=1800)
    ok = r.returncode == 0
    print(f'    validate.py rc={r.returncode}  {"✅ 通过" if ok else "❌ 失败"}')
    for l in (r.stdout or '').strip().splitlines()[-3:]:
        print('      ', l[:110])
    if not ok:
        problems.append('包内 validate.py 运行失败')
finally:
    shutil.rmtree(tmp, ignore_errors=True)

# ---------- 4. 与仓库一致性 ----------
print('\n[4] 关键文件与仓库一致性（行尾归一化后哈希）')
sample = ['docs/manuscript_v34_nar.md', 'results/benchmark/v33/permutation_v33.json',
          'scripts/rerun_v33.py', 'scripts/validate.py', 'cds_tool/cds.py',
          'requirements.txt', 'LICENSE']
diff = 0
for s in sample:
    p = REPO / s
    if not p.exists():
        notes.append(f'仓库内不存在: {s}')
        continue
    a = hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()[:12]
    b = hashlib.sha256(z.read('SPATBench/' + s).replace(b'\r\n', b'\n')).hexdigest()[:12]
    if a != b:
        diff += 1
        problems.append(f'内容不一致: {s}')
        print(f'    ❌ {s}  {a} vs {b}')
print(f'    一致 {len(sample) - diff}/{len(sample)}' + ('  ✅' if diff == 0 else ''))

# ---------- 5. 路径与编码 ----------
print('\n[5] 路径与编码')
bad_path = [n for n in names if '\\' in inner(n)]
bad_enc = [n for n in names if re.search(r'[^\x00-\x7F]', n)]
print(f'    反斜杠路径 {len(bad_path)}  {"✅" if not bad_path else "❌"}')
print(f'    非 ASCII 名 {len(bad_enc)}  {"✅" if not bad_enc else ""}')
if bad_path:
    problems.append('路径含反斜杠')

print('\n' + '=' * 64)
if problems:
    print(f'发现 {len(problems)} 个问题:')
    for p in problems:
        print('   -', p)
    sys.exit(1)
print('✅ 深度审计通过：全量格式正常、必备文件齐全、解包可运行、内容一致、路径规范')
if notes:
    for n in notes:
        print('   备注:', n)
