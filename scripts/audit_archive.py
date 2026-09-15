"""Archive integrity audit (basic).

Usage:
    python scripts/audit_archive.py [path/to/SPATBench_v33_v34_zenodo.zip]
If no argument is given, looks for SPATBench_v33_v34_zenodo.zip in the
repository root (the packaging script writes it there).

Checks: zip CRC integrity, file inventory, forbidden content, zero-byte files,
and readability sampling of JSON / PNG members.
"""
import io
import json
import random
import sys
import zipfile
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ZIP = (Path(sys.argv[1]) if len(sys.argv) > 1
       else REPO / 'SPATBench_v33_v34_zenodo.zip')

if not ZIP.exists():
    print(f'找不到归档文件: {ZIP}')
    print('用法: python scripts/audit_archive.py [path/to/SPATBench_v33_v34_zenodo.zip]')
    print('（无参数时默认在仓库根目录查找）')
    sys.exit(2)

print(f'检查目标: {ZIP.name}  {ZIP.stat().st_size / 1024 / 1024:.1f} MB\n')
problems = []

with zipfile.ZipFile(ZIP) as z:
    bad = z.testzip()
    if bad:
        problems.append(f'CRC 校验失败: {bad}')
        print(f'[1] CRC 校验: 损坏文件 {bad}')
    else:
        print('[1] CRC 校验: 通过（无损坏）')
    names = z.namelist()
    infos = [z.getinfo(n) for n in names if not n.endswith('/')]
    print(f'    文件数: {len(infos)}')

    tops = Counter()
    for i in infos:
        parts = i.filename.split('/')
        tops['/'.join(parts[:2]) if len(parts) > 2 else parts[0]] += 1
    print('\n[2] 目录构成:')
    for k, v in tops.most_common(15):
        print(f'    {v:4d}  {k}')

    forbid = {
        '零字节文件': [i.filename for i in infos if i.file_size == 0],
        '历史草稿 FINAL': [i.filename for i in infos if '_FINAL' in i.filename],
        '调试图': [i.filename for i in infos if any(
            k in i.filename for k in ('fitz_', 'checkbox_zoom', 'titlepage_',
                                      'di_form_p', 'di_visible', 'check_verify'))],
        '补丁脚本': [i.filename for i in infos if 'patch_manuscript' in i.filename],
        '虚拟环境/缓存': [i.filename for i in infos if any(
            k in i.filename for k in ('venv/', '__pycache__', '.egg-info',
                                      '.pytest_cache'))],
        '过程日志': [i.filename for i in infos if i.filename.endswith('.log')],
        # 内部评审/交接/运维材料 —— 永不属于公开交付集（与 pack_zenodo 的
        # INTERNAL_PATHS、adversarial_audit 的 EXCL_INTERNAL 保持一致）
        '内部评审与运维材料': [i.filename for i in infos if any(
            k in i.filename for k in (
                'docs/reviews/', 'HANDOVER_2026-09-04', 'zenodo_upload_guide',
                'github_release_guide', 'data_repository_alternatives',
                'declaration_of_interests_guide', 'adversarial_audit_report',
                'deep_audit_report', 'script_verification_log', 'watchdog_gide',
                'run_gide_extend', 'push_to_github', 'prepare_github_repos',
                'normalise_repos', 'rewrite_github_history',
                'sync_github_deterministic', 'README_v33_summary',
                'run_all_checks_report', 'cover_letter',
                'di_form_filled'))],
    }
    print('\n[3] 禁止内容检查:')
    for label, hits in forbid.items():
        if hits:
            problems.append(f'{label}: {len(hits)} 个')
            print(f'    ❌ {label}: {len(hits)} 个  {hits[:2]}')
        else:
            print(f'    ✅ {label}: 0')

    key = {
        '数据 h5ad': [i.filename for i in infos if i.filename.endswith('.h5ad')],
        '结果 JSON': [i.filename for i in infos if i.filename.endswith('.json')],
        '图 PNG': [i.filename for i in infos if i.filename.endswith('.png')],
        '脚本 py': [i.filename for i in infos if i.filename.endswith('.py')],
    }
    print('\n[4] 关键内容:')
    for k, v in key.items():
        print(f'    {k}: {len(v)}')

    print('\n[5] 抽样可读性:')
    jsons = [i.filename for i in infos
             if i.filename.endswith('.json') and i.file_size > 200]
    for n in random.sample(jsons, min(4, len(jsons))):
        try:
            json.loads(z.read(n))
            print(f'    ✅ JSON: {n.split("/", 1)[-1][:60]}')
        except Exception as e:
            problems.append(f'JSON 损坏: {n}')
            print(f'    ❌ JSON: {n} -> {e}')
    pngs = [i.filename for i in infos if i.filename.endswith('.png')]
    for n in random.sample(pngs, min(2, len(pngs))):
        try:
            from PIL import Image
            im = Image.open(io.BytesIO(z.read(n)))
            im.load()
            print(f'    ✅ PNG {im.size[0]}x{im.size[1]}: {n.split("/", 1)[-1][:50]}')
        except Exception as e:
            problems.append(f'PNG 损坏: {n}')
            print(f'    ❌ PNG: {n} -> {e}')

print('\n' + '=' * 60)
if problems:
    print(f'发现 {len(problems)} 个问题:')
    for p in problems:
        print('   -', p)
    sys.exit(1)
print('✅ 全部通过：CRC 完整、无禁止内容、抽样文件可读')
