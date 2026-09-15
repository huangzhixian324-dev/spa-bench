"""Multiple-comparison family sensitivity analysis (review P0 item 1).

Recomputes BH q-values for the complete 36-cell permutation matrix under
three family definitions and writes a machine-readable JSON plus a
markdown table for Supplementary Table S21:

  A. within-cohort (manuscript primary): 6 families of 6 cells each
  B. global single family: 1 family of 36 cells
  C. stratified by method class: fixed scorers (24 cells) and
     trainable (12 cells)

Purpose: show explicitly which significance calls depend on the family
choice (reviewer concern R2-M3), so the manuscript can disclose it.
"""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
V33 = REPO / "results" / "benchmark" / "v33"


def bh(pairs):
    pairs = sorted(pairs, key=lambda t: t[1])
    m = len(pairs)
    prev, out = 1.0, {}
    for rank in range(m, 0, -1):
        k, p = pairs[rank - 1]
        prev = min(prev, p * m / rank, 1.0)
        out[k] = prev
    return out


def main():
    perms = json.load(open(V33 / "permutation_v33.json"))
    cells = {}
    for c, methods in perms["cohorts"].items():
        for m, d in methods.items():
            cells[(c, m)] = d["p"]
    assert len(cells) == 36, f"expected 36 cells, got {len(cells)}"

    fixed = {"IMPRES", "GEP", "TIDE", "PD_L1"}
    fams = {
        "A_within_cohort": {},
        "B_global": bh(list(cells.items())),
        "C_stratified": {},
    }
    for c in perms["cohorts"]:
        fam = {k: p for k, p in cells.items() if k[0] == c}
        fams["A_within_cohort"][c] = bh(list(fam.items()))
    fams["C_stratified"]["fixed_24"] = bh(list(
        {k: p for k, p in cells.items() if k[1] in fixed}.items()))
    fams["C_stratified"]["trainable_12"] = bh(list(
        {k: p for k, p in cells.items() if k[1] not in fixed}.items()))

    rows = []
    for (c, m), p in sorted(cells.items(), key=lambda t: t[1]):
        qA = fams["A_within_cohort"][c][(c, m)]
        qB = fams["B_global"][(c, m)]
        qC = (fams["C_stratified"]["fixed_24"].get((c, m))
              if m in fixed
              else fams["C_stratified"]["trainable_12"].get((c, m)))
        sig = lambda q: "Yes" if q < 0.05 else ""
        rows.append({"cohort": c, "method": m, "p": p,
                     "q_within": qA, "q_global": qB, "q_strat": qC,
                     "sig_within": sig(qA), "sig_global": sig(qB),
                     "sig_strat": sig(qC)})

    out = {"model": "BH adjusted p-values over the 36-cell permutation "
                    "matrix under three family definitions",
           "families": {"A": "within-cohort, 6 families x 6 cells "
                             "(manuscript primary)",
                        "B": "global, 1 family x 36 cells",
                        "C": "stratified by method class: fixed scorers "
                             "24 cells / trainable 12 cells"},
           "rows": rows}
    with open(V33 / "bh_sensitivity_v33.json", "w") as fh:
        json.dump(out, fh, indent=1)

    # summary of calls that change
    changed = [r for r in rows if r["sig_within"] != r["sig_global"]
               or r["sig_within"] != r["sig_strat"]]
    print(f"{'cohort':22s} {'method':14s} {'p':>9s} "
          f"{'qA':>7s} {'qB':>7s} {'qC':>7s}  sig(A/B/C)")
    for r in rows:
        mark = "  <-- differs" if r in changed else ""
        print(f"{r['cohort']:22s} {r['method']:14s} {r['p']:9.5f} "
              f"{r['q_within']:7.4f} {r['q_global']:7.4f} "
              f"{r['q_strat']:7.4f}   "
              f"{r['sig_within'] or '-'}/{r['sig_global'] or '-'}/"
              f"{r['sig_strat'] or '-'}{mark}")
    print(f"\n{len(changed)} cells change significance across families")
    print("written", V33 / "bh_sensitivity_v33.json")


if __name__ == "__main__":
    main()
