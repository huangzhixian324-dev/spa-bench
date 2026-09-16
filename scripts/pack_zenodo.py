"""Package the complete SPATBench artifact for Zenodo upload.
Includes code, data, results, docs, CDS package; excludes .git / venv / caches
and all internal process material (INTERNAL_PATHS).
"""
import zipfile
from pathlib import Path

R = Path(__file__).resolve().parents[1]
REPO = R  # alias
OUT = REPO / 'SPATBench_v33_v34_zenodo.zip'
EXCLUDE_DIRS = {'.git', 'venv', '__pycache__', 'github-upload', '.pytest_cache',
                'node_modules', '_build', 'archive', '_process'}
EXCLUDE_FILES = {'SPATBench_v33_v34_zenodo.zip'}
EXCLUDE_SUFFIX = {'.log', '.pyc'}
EXCLUDE_NAME_PARTS = ('egg-info',)

# Internal process material — must NEVER be part of the public deliverable:
# peer-review simulations, development handover notes, local upload how-tos,
# internal audit/verification reports and personal machine tooling. They live
# in the working tree or in the local archive only. The same list is mirrored
# in scripts/adversarial_audit.py and in the local (unpublished) copy of
# sync_github_deterministic.py.
INTERNAL_PATHS = (
    'docs/reviews',
    'docs/HANDOVER_2026-09-04.md',
    'docs/zenodo_upload_guide.md',
    'docs/github_release_guide.md',
    'docs/data_repository_alternatives.md',
    'docs/declaration_of_interests_guide.md',
    # Submission correspondence — editor-confidential / personal data, never public:
    'docs/cover_letter_cellsystems.md',
    'docs/di_form_filled.pdf',
    'results/adversarial_audit_report.txt',
    'results/deep_audit_report.txt',
    'results/script_verification_log.txt',
    'scripts/watchdog_gide.py', 'scripts/watchdog_gide.ps1', 'scripts/watchdog_gide.sh',
    'scripts/run_gide_extend.bat', 'scripts/push_to_github.py',
    'scripts/prepare_github_repos.py', 'scripts/normalise_repos.py',
    'scripts/rewrite_github_history.py', 'scripts/sync_github_deterministic.py',
    'results/benchmark/v33/README_v33_summary.md',
    'results/run_all_checks_report.txt',
)


def _internal(rel: str) -> bool:
    return any(rel == p or rel.startswith(p + '/') for p in INTERNAL_PATHS)


n = 0
total = 0
with zipfile.ZipFile(OUT, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for f in sorted(R.rglob('*')):
        if not f.is_file():
            continue
        rel = f.relative_to(R).as_posix()
        if _internal(rel):
            continue
        parts = set(f.relative_to(R).parts)
        if parts & EXCLUDE_DIRS:
            continue
        if f.name in EXCLUDE_FILES:
            continue
        if f.suffix in EXCLUDE_SUFFIX:
            continue
        if any(p in f.name for p in EXCLUDE_NAME_PARTS):
            continue
        z.write(f, f'SPATBench/{rel}')
        n += 1
        total += f.stat().st_size
print(f'打包完成: {n} 文件, 原始 {total/1024/1024:.1f} MB -> zip {OUT.stat().st_size/1024/1024:.1f} MB')
print('位置:', OUT)
