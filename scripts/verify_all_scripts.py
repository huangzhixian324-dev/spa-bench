"""Run and verify every remaining runnable script; log results.
Tier A: standalone result/deliverable scripts (run for real).
Tier B: data-dependent scripts (compile-verified already; run those that are quick).
Long permutation jobs (rerun_v33, fill_nc_cells, run_null_200) are reported as
intentionally not re-run (already completed and archived).
"""
import subprocess
import time
from pathlib import Path

R = Path(__file__).resolve().parents[1]
PY = str(R) + '/venv/Scripts/python.exe'
LOG = R / 'results' / 'script_verification_log.txt'

TIER_A = [
    ('validate.py', 600),
    ('power_table_v33.py', 300),
    ('bh_sensitivity_v33.py', 300),
    ('mi_var_overlap_v33.py', 300),
    ('hugo_envar_clinical_v33.py', 300),
    ('make_run_registry.py', 300),
    ('integrity_audit.py', 600),
    ('make_figures_v33.py', 900),
    ('make_graphical_abstract.py', 300),
    ('synthesize_v33_supplementary.py', 600),
    ('build_submission.py', 300),
    ('md_to_docx_cellpress.py', 300),
    ('pack_zenodo.py', 900),
]
TIER_B_QUICK = [
    ('cds_negative_controls.py', 600),
    ('e3_literature_backtest.py', 600),
    ('clinical_utility_and_synth_v33.py', 600),
]
SKIP_LONG = {
    'rerun_v33.py': '全量置换管线（已跑完 3,200+ 次，结果已存档）',
    'fill_nc_cells.py': '置换扩展（Gide 5,000 次已跑完并归档）',
    'run_null_200.py': '200 复本 null（已完成）',
    'extend_final.py': 'Gide 5,000 次续跑（已完成）',
    'bench_gide.py': '旧 Gide 基准（被 v33 取代）',
    'rerun_v33.py': '同上',
    'e1_gse274975_replication.py': '被取消资格的队列（记录在案）',
    'download_all.py': '数据下载（数据已在本地）',
    'download_all_cohorts.py': '同上',
    'download_cohorts.py': '同上',
    'download_cohorts_ena.py': '同上',
    'download_robust.py': '同上',
    'download_tcga.py': '同上',
    'fetch_cbioportal_cohorts.py': '同上',
    'fetch_gide_api.py': '同上（API 抓取）',
    'fetch_gide_full.py': '同上',
    'rebuild_gide.py': '数据重建（已完成）',
    'rebuild_gide_splits.py': '同上',
    'rebuild_jung.py': '同上',
    'rebuild_nathanson_splits.py': '同上',
    'rebuild_riaz_recist.py': '同上（记录在修正日志）',
    'map_riaz_to_hgnc.py': '数据映射（已完成）',
    'verify_restored_data.py': '数据校验（已完成）',
    'preprocess_bulk.py': '预处理（已完成）',
    'build_coad_gdc.py': 'TCGA-COAD 重获取（已完成）',
    'run_tcga_cox.py': 'TCGA Cox（已完成）',
    'build_graph.py': '图构建辅助（一次性）',
    'catalog_cohorts.py': '清点（一次性）',
    'catalog_spatial.py': '清点（一次性）',
    'generate_real_figures.py': '旧图脚本（被 make_figures_v33 取代）',
    'supplementary_analysis.py': '旧补充分析（被 v33 取代）',
    'supplementary_v33_experiments.py': '补充实验（已完成）',
    'final_benchmarks_v33.py': '基准汇总（已完成）',
    'run_full_pipeline.py': '全管线入口（含长任务）',
    'merge_nc_cells.py': 'n.c. 合并（已完成）',
    'pdl1_perm_v33.py': 'PD-L1 置换（已完成）',
    'impres_perm_riaz.py': 'IMPRES 置换（已完成）',
    'power_and_xgb_archive.py': 'one-sample power + XGBoost 归档（已完成，JSON 已入库）',
    'bootstrap_effect_ratio.py': '效应量 bootstrap（已完成，需先跑 cv_predict；JSON 已入库）',
    'rebuild_refs_cellpress2.py': '引用重建（已完成，稿件已更新）',
    'watchdog_gide.py': 'Gide 看门狗（服务型脚本，长驻）',
    'e1_replication_imvigor210.py': '复制实验（已完成）',
    'e1_replication_liu2019.py': '复制实验（已完成）',
    'hugo_survival_true.py': '生存分析（已完成）',
    'hugo_survival_v33.py': '生存分析（已完成）',
    'cds_null_all_endpoints.py': 'CDS null（已完成，200 复本）',
}

lines = []
def log(s):
    print(s, flush=True)
    lines.append(s)

log(f'=== 脚本验证 {time.strftime("%Y-%m-%d %H:%M")} ===\n')
log('--- Tier A（实跑，交付物/结果脚本）---')
for script, tmo in TIER_A:
    p = R / 'scripts' / script
    t0 = time.time()
    try:
        r = subprocess.run([PY, str(p)], cwd=str(R), capture_output=True,
                           text=True, timeout=tmo)
        dt = time.time() - t0
        err = (r.stderr or '').strip().splitlines()
        log(f'  [{"OK" if r.returncode==0 else "RC="+str(r.returncode)}] {script}  {dt:.1f}s'
            + (f'  | {err[-1][:80]}' if err and r.returncode != 0 else ''))
    except subprocess.TimeoutExpired:
        log(f'  [TIMEOUT] {script} >{tmo}s')

log('\n--- Tier B（快速数据脚本）---')
for script, tmo in TIER_B_QUICK:
    p = R / 'scripts' / script
    t0 = time.time()
    try:
        r = subprocess.run([PY, str(p)], cwd=str(R), capture_output=True,
                           text=True, timeout=tmo)
        dt = time.time() - t0
        err = (r.stderr or '').strip().splitlines()
        log(f'  [{"OK" if r.returncode==0 else "RC="+str(r.returncode)}] {script}  {dt:.1f}s'
            + (f'  | {err[-1][:80]}' if err and r.returncode != 0 else ''))
    except subprocess.TimeoutExpired:
        log(f'  [TIMEOUT] {script} >{tmo}s')

log('\n--- 未重跑（已完成/长驻/需网络，附理由）---')
for name, why in sorted(SKIP_LONG.items()):
    if not (R / 'scripts' / name).exists():
        continue          # moved to the process archive; no longer part of the tree
    log(f'  {name}: {why}')
archived = sorted(n for n in SKIP_LONG if not (R / 'scripts' / n).exists())
if archived:
    log(f'\n（另有 {len(archived)} 个脚本已移出交付集、存入过程存档，不在本清单内）')

LOG.write_text('\n'.join(lines), encoding='utf-8')
print('\nlog saved:', LOG)
