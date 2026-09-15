"""One-command verification of the SPATBench archive.

Run from the repository root (or anywhere):
    python scripts/run_all_checks.py

Checks, in order:
  1. dependency imports
  2. end-to-end self-test        (scripts/validate.py, synthetic cohort)
  3. power analysis              (scripts/power_table_v33.py)
  4. FDR family sensitivity      (scripts/bh_sensitivity_v33.py)
  5. feature-selection overlap   (scripts/mi_var_overlap_v33.py)
  6. Hugo EN-Var clinical cell   (scripts/hugo_envar_clinical_v33.py)
  7. figure regeneration         (scripts/make_figures_v33.py)
  8. graphical abstract          (scripts/make_graphical_abstract.py)
  9. supplementary regeneration  (scripts/synthesize_v33_supplementary.py)
 10. submission PDF assembly     (scripts/build_submission.py)
 11. Word main text             (scripts/md_to_docx_cellpress.py)
 12. archive integrity audit     (scripts/audit_archive.py)

The script locates the repository from its own path, so it works after
unzipping the Zenodo archive anywhere. Exit code 0 = everything passed.
"""
import importlib
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
REPO = HERE.parents[1]
PY = sys.executable

DEPENDENCIES = ['numpy', 'scipy', 'pandas', 'sklearn', 'anndata', 'h5py',
                'tidepy', 'lifelines', 'matplotlib', 'PIL', 'markdown',
                'docx', 'pypdf']

SCRIPTS = [
    ('端到端自检（合成数据全管线）', 'validate.py', 1200),
    ('功效分析', 'power_table_v33.py', 300),
    ('FDR 家族敏感性', 'bh_sensitivity_v33.py', 300),
    ('特征选择重叠', 'mi_var_overlap_v33.py', 600),
    ('Hugo EN-Var 临床格', 'hugo_envar_clinical_v33.py', 600),
    ('全部图（300 dpi）', 'make_figures_v33.py', 900),
    ('图形摘要', 'make_graphical_abstract.py', 300),
    ('补充材料再生', 'synthesize_v33_supplementary.py', 600),
    ('投稿 PDF 组装', 'build_submission.py', 300),
    ('Word 主稿', 'md_to_docx_cellpress.py', 300),
    ('存档完整性审计', 'audit_archive.py', 900),
]


def main() -> int:
    print('=' * 70)
    print(f'SPATBench 验证  |  仓库: {REPO}')
    print(f'Python: {sys.version.split()[0]}  |  {PY}')
    print('=' * 70)
    results = []

    # 1) dependencies
    print('\n[1/12] 依赖检查')
    missing = []
    for m in DEPENDENCIES:
        try:
            importlib.import_module(m)
            print(f'    ✅ {m}')
        except Exception as e:
            missing.append(m)
            print(f'    ❌ {m}  ({type(e).__name__})')
    results.append(('依赖检查', 'PASS' if not missing else 'FAIL',
                    0.0, f'缺失: {missing}' if missing else '全部可用'))
    if missing:
        print(f'\n缺少依赖: {missing}\n请先运行: pip install -r requirements.txt')
        return 1

    # 2-12) scripts
    for idx, (label, script, tmo) in enumerate(SCRIPTS, start=2):
        p = REPO / 'scripts' / script
        print(f'\n[{idx}/12] {label}  ({script})')
        if not p.exists():
            print('    ❌ 脚本不存在')
            results.append((label, 'FAIL', 0.0, '脚本缺失'))
            continue
        t0 = time.time()
        try:
            r = subprocess.run([PY, str(p)], cwd=str(REPO), capture_output=True,
                               text=True, timeout=tmo)
            dt = time.time() - t0
            ok = r.returncode == 0
            out = (r.stdout or '').strip().splitlines()
            tail = out[-1][:90] if out else ''
            print(f'    {"✅ PASS" if ok else "❌ FAIL"}  {dt:.1f}s  {tail}')
            if not ok:
                err = (r.stderr or '').strip().splitlines()
                print('    stderr:', (err[-1][:150] if err else '(none)'))
            results.append((label, 'PASS' if ok else f'FAIL(rc={r.returncode})',
                            dt, tail))
        except subprocess.TimeoutExpired:
            print(f'    ⏱ TIMEOUT (>{tmo}s)')
            results.append((label, 'TIMEOUT', float(tmo), ''))

    # summary
    print('\n' + '=' * 70)
    print('汇总')
    print('=' * 70)
    npass = sum(1 for _, s, _, _ in results if s == 'PASS')
    for label, status, dt, note in results:
        print(f'  {status:12s} {dt:7.1f}s  {label}  {note[:50]}')
    print(f'\n通过 {npass}/{len(results)}')
    ok_all = npass == len(results)
    print('✅ 全部通过——存档可复现' if ok_all else '⚠️ 存在未通过项（见上）')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
