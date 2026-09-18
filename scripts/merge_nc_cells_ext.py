"""Merge the 8 extended permutation cells (dry experiment A) into the
permutation matrix, following the established merge_nc_cells.py convention.

- backs up permutation_v33.json first;
- overlays p / n_shuffle / obs_auroc for the 8 n<=43 trainable cells from
  results/benchmark/v33/nc_cells_ext/*.json (each records its archived
  values for the audit trail);
- copies the ext JSONs into nc_cells/ so the supplementary generator's
  Perm-obs lookup picks up the new frozen-protocol observations;
- prints the frozen-vs-benchmark deviations for the S11 footnote.
"""
import json
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
V33 = REPO / "results" / "benchmark" / "v33"
PMF = V33 / "permutation_v33.json"
EXT = V33 / "nc_cells_ext"
NCD = V33 / "nc_cells"

CELLS = [
    ("Hugo_2016", "ElasticNet"), ("Hugo_2016", "ElasticNet_Var"),
    ("Nathanson_2017", "ElasticNet"), ("Nathanson_2017", "ElasticNet_Var"),
    ("Jung_2019_DCB", "ElasticNet"), ("Jung_2019_DCB", "ElasticNet_Var"),
    ("Riaz_2017_RECIST_v33", "ElasticNet"),
    ("Riaz_2017_RECIST_v33", "ElasticNet_Var"),
]


def main():
    bak = V33 / "permutation_v33.pre_ext_backup.json"
    shutil.copy2(PMF, bak)
    pm = json.load(open(PMF))
    BM = json.load(open(V33 / "benchmark_v33.json"))["cohorts"]

    devs = []
    for cohort, method in CELLS:
        ext = json.load(open(EXT / f"{cohort}__{method}.json"))
        cell = pm["cohorts"][cohort][method]
        old = {"p": cell.get("p"), "n_shuffle": cell.get("n_shuffle"),
               "obs_auroc": cell.get("obs_auroc")}
        cell["p"] = ext["p"]
        cell["n_shuffle"] = ext["n_perm"]
        cell["obs_auroc"] = ext["obs_auroc_this_run"]
        cell["ext_2026_09_18"] = {
            "archived": old, "n_valid": ext["n_valid"],
            "protocol": ext["protocol"]}
        bm_a = BM[cohort][method]["auroc"]
        devs.append((f"{cohort}/{method}", ext["obs_auroc_this_run"],
                     bm_a, abs(ext["obs_auroc_this_run"] - bm_a)))
        # sync nc_cells so the generator's Perm-obs lookup sees the new obs
        shutil.copy2(EXT / f"{cohort}__{method}.json",
                     NCD / f"{cohort}__{method}.json")

    tmp = str(PMF) + ".tmp"
    json.dump(pm, open(tmp, "w"), indent=1)
    shutil.move(tmp, PMF)

    print("merged 8 extended cells; deviations (frozen obs vs benchmark):")
    for name, obs, bm, d in sorted(devs, key=lambda t: -t[3]):
        print(f"  {name:38} obs={obs:.4f} bm={bm:.4f} dev={d:.4f}")
    print("backup:", bak)
    return 0


if __name__ == "__main__":
    sys.exit(main())
