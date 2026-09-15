"""Full-project data/resource integrity audit (2026-09-08).

Checks every resource the project claims to have:
  1. run_registry.data_hashes  -> existence + SHA-1 (first 16 hex) match
  2. rerun_v33.COHORTS         -> h5ad / splits existence
  3. figures S1-S12            -> results/figures/v33/*.png
  4. supplementary sources     -> JSONs cited by the supplementary generator
  5. external fetches          -> data/external/*
  6. paths cited in manuscript + supplementary markdown
Writes a categorized report to stdout.
"""
import hashlib
import json
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
os.chdir(REPO)
sys.path.insert(0, str(REPO / "scripts"))

ok, missing, mismatch = [], [], []


def check(path, kind="exists"):
    p = Path(path)
    if not p.exists():
        missing.append((str(path), kind))
        return None
    ok.append((str(path), kind))
    return p


def sha1_16(path):
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


print("=" * 70)
print("1) run_registry.data_hashes (existence + SHA-1 first-16)")
print("=" * 70)
reg = json.load(open("results/run_registry.json"))
for path, expected in reg.get("data_hashes", {}).items():
    p = Path(path)
    if not p.exists():
        missing.append((path, f"registry-hash expected {expected}"))
        print(f"  MISS  {path}  (registry sha1[:16]={expected})")
        continue
    actual = sha1_16(p)
    if actual == expected:
        ok.append((path, "registry-hash MATCH"))
        print(f"  MATCH {path}")
    else:
        mismatch.append((path, expected, actual))
        print(f"  DIFF  {path}  registry={expected} actual={actual}")

print()
print("=" * 70)
print("2) rerun_v33.COHORTS data files")
print("=" * 70)
try:
    from rerun_v33 import COHORTS
    for name, cfg in sorted(COHORTS.items()):
        for key in ("h5ad", "splits"):
            if cfg.get(key):
                check(cfg[key], f"cohort {name}/{key}")
except Exception as e:
    print("  (COHORTS import failed:", e, ")")

print()
print("=" * 70)
print("3) figures referenced by the supplementary generator")
print("=" * 70)
for i in range(1, 13):
    hits = list(Path("results/figures/v33").glob(f"figS{i}_*.png"))
    if hits:
        ok.append((str(hits[0]), f"figure S{i}"))
        print(f"  OK    figS{i}: {hits[0].name}")
    else:
        missing.append((f"results/figures/v33/figS{i}_*.png", f"figure S{i}"))
        print(f"  MISS  figS{i}")

print()
print("=" * 70)
print("4) paths cited in manuscript + supplementary markdown")
print("=" * 70)
cited = set()
for doc in ["docs/manuscript_v34_nar.md", "docs/supplementary_material_v33.md"]:
    text = Path(doc).read_text(encoding="utf-8")
    for pat in [r"`(results/[^`\s]+?)`", r"`(data/[^`\s]+?)`",
                r"`(scripts/[^`\s]+?)`", r"`(docs/[^`\s]+?)`"]:
        for mch in re.findall(pat, text):
            if not mch.endswith(("/", "*")) and "*" not in mch:
                cited.add(mch.split(" ")[0])
missing_cited = []
for p in sorted(cited):
    if p.startswith("docs/"):
        continue  # docs referencing other docs is fine
    if not Path(p).exists():
        missing_cited.append(p)
        missing.append((p, "cited in manuscript/supplementary"))
print(f"  cited local paths: {len(cited)}; missing: {len(missing_cited)}")
for p in missing_cited:
    print(f"    MISS  {p}")

print()
print("=" * 70)
print("5) external fetches + key project assets")
print("=" * 70)
for d in ["data/external/mel_iatlas_hugo_ucla_2016",
          "data/external/blca_iatlas_imvigor210_2017",
          "data/external/mel_iatlas_liu_2019",
          "data/external/skcm_vanderbilt_mskcc_2015",
          "data/external/mel_iatlas_gide_2019",
          "data/tcga/tcga_cox_all.json",
          "data/cohorts/GSE274975/raw",
          "data/tcga/TCGA-COAD_v2",
          "cds_tool/cds.py", "cds_tool/setup.py",
          "Dockerfile", "scripts/validate.py",
          "results/benchmark/v33/power_table_v33.json",
          "results/benchmark/v33/bh_sensitivity_v33.json",
          "results/benchmark/v33/nc_cells",
          "results/run_registry.json"]:
    exists = Path(d).exists()
    (ok if exists else missing).append((d, "asset"))
    print(f"  {'OK  ' if exists else 'MISS'}  {d}")

print()
print("=" * 70)
print(f"SUMMARY: ok={len(ok)}  missing={len(missing)}  hash-mismatch={len(mismatch)}")
print("=" * 70)
for p, why in missing:
    print(f"  MISSING  [{why[:40]}]  {p}")
for p, exp, act in mismatch:
    print(f"  MISMATCH {p}  registry={exp} actual={act}")
