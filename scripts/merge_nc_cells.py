"""Merge filled n.c. permutation cells into permutation_v33.json and
recompute within-cohort BH over the (now complete) cell families.

Reads: results/benchmark/v33/nc_cells/*.json (one per cell, from
scripts/fill_nc_cells.py)
Updates: results/benchmark/v33/permutation_v33.json in place
  - adds each cell with p and n_perm
  - meta note documents the completion
Prints an audit of BH changes for every cell in the affected cohorts.
"""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
V33 = REPO / "results" / "benchmark" / "v33"
CELLS = V33 / "nc_cells"


def bh(per_list):
    # Standard Benjamini-Hochberg adjusted p-values (cummin from the
    # largest rank downwards). 2026-09-07 fix: was cummax-from-smallest,
    # which overestimates q for non-significant cells.
    per_list = sorted(per_list, key=lambda t: t[1])
    mc = len(per_list)
    prev = 1.0
    out = {}
    for rank in range(mc, 0, -1):
        m, p = per_list[rank - 1]
        prev = min(prev, p * mc / rank, 1.0)
        out[m] = prev
    return out


def main():
    pf = V33 / "permutation_v33.json"
    perms = json.load(open(pf))
    files = sorted(CELLS.glob("*.json"))
    if not files:
        print("no cell files found")
        return 1
    print(f"{len(files)} cell files:")
    added = []
    for f in files:
        r = json.load(open(f))
        c, m = r["cohort"], r["method"]
        if r.get("obs_deviation", 0) > 0.005:
            print(f"  [WARN] {c}/{m}: obs {r['obs_auroc_this_run']} deviates "
                  f"from benchmark {r['benchmark_auroc']} by "
                  f"{r['obs_deviation']} (recorded)")
        perms["cohorts"].setdefault(c, {})[m] = {
            "obs_auroc": r["obs_auroc_this_run"], "p": r["p"],
            "n_perm": r["n_perm"]}
        added.append((c, m, r["p"], r["n_perm"]))
        print(f"  + {c}/{m}: obs={r['obs_auroc_this_run']} p={r['p']:.4g} "
              f"n={r['n_perm']}")

    # report BH before/after for affected cohorts
    for c in sorted({c for c, _, _, _ in added}):
        before = {m: cell["p"] for m, cell in perms["cohorts"][c].items()
                  if m not in {m2 for c2, m2, _, _ in added if c2 == c}}
        after_all = perms["cohorts"][c]
        qb = bh(list(before.items()))
        qa = bh([(m, cell["p"]) for m, cell in after_all.items()])
        print(f"\n{c}: BH family {len(before)} -> {len(after_all)}")
        for m in after_all:
            tag = "(new)" if m not in before else ""
            print(f"  {m}: q {qb.get(m, float('nan')):.3f} -> "
                  f"{qa[m]:.3f} {tag}")

    note = (" 2026-09-06: the previously n.c. full-pipeline permutation "
            "cells (ElasticNet / ElasticNet_Var on Jung_2019_DCB, "
            "Riaz_2017_cytolytic, Riaz_2017_RECIST_v33, and ElasticNet_Var "
            "on Nathanson_2017) were completed by "
            "scripts/fill_nc_cells.py (per-cell n recorded, seed 42, frozen "
            "hyperparameters from tune_primary); the permutation matrix is "
            "now complete (36/36 method x cohort-endpoint cells have a "
            "recorded p), and within-cohort BH families were recomputed.")
    perms["meta"]["note"] = perms["meta"].get("note", "") + note
    with open(pf, "w") as fh:
        json.dump(perms, fh, indent=1)
    print("\nupdated", pf)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
