"""Hugo 2016 survival with TRUE OS events (replaces the OS<=400d approximation).

Requires the iAtlas harmonization of Hugo 2016 fetched by
scripts/e3_literature_backtest.py (data/external/mel_iatlas_hugo_ucla_2016/).

Pipeline:
  1. match the 28 local Hugo rows (no barcodes) to iAtlas samples by
     per-sample expression-profile correlation with a GREEDY one-to-one
     assignment (best pair first), reporting margins;
  2. cross-check: local OS_days vs iAtlas OS_MONTHS x 30.436875
     (agreement corroborates the matching);
  3. re-run the survival analysis with true events
     (event = OS_STATUS '1:DECEASED', duration = OS_MONTHS):
     IMPRES (canonical) and GEP z-scored Cox per-SD HR, median-split KM
     log-rank, RECIST KM log-rank.

Output: results/benchmark/v33/hugo_survival_true.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
EXT = REPO / "data" / "external" / "mel_iatlas_hugo_ucla_2016"
RES = REPO / "results" / "benchmark" / "v33"
DAYS_PER_MONTH = 30.436875
GEP_GENES = ["CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1",
             "TBX21"]


def main():
    mat = pd.read_csv(EXT / "expression.tsv.gz", index_col=0)
    clin = pd.read_csv(EXT / "clinical.tsv", sep="\t", index_col=0)
    import scanpy as sc
    a = sc.read_h5ad(REPO / "data" / "cohorts" / "Hugo_2016" / "processed" /
                     "Hugo_2016_processed.h5ad")
    X = (a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X))
    genes = list(a.var_names)
    common = [g for g in mat.index if g in genes]
    li = {g: i for i, g in enumerate(genes)}
    A = X[:, [li[g] for g in common]]
    B = mat.loc[common].T.values
    Az = (A - A.mean(1, keepdims=True)) / A.std(1, keepdims=True)
    Bz = (B - B.mean(1, keepdims=True)) / B.std(1, keepdims=True)
    corr = Az @ Bz.T / len(common)

    # greedy 1-1 assignment
    n_loc, n_ext = corr.shape
    pairs = sorted(((corr[i, j], i, j) for i in range(n_loc)
                    for j in range(n_ext)), reverse=True)
    used_i, used_j, assign = set(), set(), {}
    for c, i, j in pairs:
        if i in used_i or j in used_j:
            continue
        used_i.add(i)
        used_j.add(j)
        assign[i] = (j, c)
        if len(assign) == min(n_loc, n_ext):
            break

    rows = []
    for i in range(n_loc):
        if i not in assign:
            rows.append({"local_row": i, "iatlas_sample": None,
                         "corr": None, "margin": None,
                         "os_status": None, "os_months": None,
                         "response_iatlas": None,
                         "local_os_days": float(a.obs["OS_days"].iloc[i]),
                         "local_response": int(a.obs["response"].iloc[i]),
                         "note": "no iAtlas mRNA sample available "
                                 "(one sample lacks expression)"})
            continue
        j, c = assign[i]
        pid = mat.columns[j]
        second = max(corr[i, jj] for jj in range(n_ext) if jj != j)
        row = clin.loc[pid] if pid in clin.index else {}
        rows.append({
            "local_row": i, "iatlas_sample": pid,
            "corr": round(float(c), 3),
            "margin": round(float(c - second), 3),
            "os_status": row.get("OS_STATUS"),
            "os_months": (None if row.get("OS_MONTHS") in (None, "", "NA")
                          else float(row.get("OS_MONTHS"))),
            "response_iatlas": row.get("RESPONSE"),
            "local_os_days": float(a.obs["OS_days"].iloc[i]),
            "local_response": int(a.obs["response"].iloc[i]),
        })
    agree = [r for r in rows if r["os_months"] is not None and
             abs(r["local_os_days"] - r["os_months"] * DAYS_PER_MONTH) < 60]
    n_os = sum(1 for r in rows if r["os_months"] is not None)
    events = sum(1 for r in rows if r["os_status"] and
                 "DECEASED" in str(r["os_status"]).upper())

    # ---- survival with true events ----
    from lifelines import CoxPHFitter
    from lifelines.statistics import logrank_test
    keep = [r for r in rows if r["os_months"] is not None and
            r["os_months"] > 0]
    dur = np.array([r["os_months"] for r in keep])
    ev = np.array([1 if "DECEASED" in str(r["os_status"]).upper() else 0
                   for r in keep])
    Xk = X[[r["local_row"] for r in keep]]

    def impres_score(Xm):
        pairs_ = [("CD274", "C10orf54"), ("CD86", "CD200"), ("CD40", "CD274"),
                  ("CD28", "CD276"), ("CD40", "CD28"), ("TNFRSF14", "CD86"),
                  ("CD27", "PDCD1"), ("CD28", "CD86"), ("CD40", "CD80"),
                  ("CD40", "PDCD1"), ("CD80", "TNFSF9"), ("CD86", "HAVCR2"),
                  ("CD86", "TNFSF4"), ("CTLA4", "TNFSF4"),
                  ("PDCD1", "TNFSF4")]
        s = np.zeros(Xm.shape[0])
        m = 0
        for ga, gb in pairs_:
            if ga in genes and gb in genes:
                s += (Xm[:, genes.index(ga)] > Xm[:, genes.index(gb)])
                m += 1
        return s * (15.0 / max(m, 1)) / 15.0

    def gep_score(Xm):
        idx = [genes.index(g) for g in GEP_GENES if g in genes]
        return Xm[:, idx].mean(axis=1)

    out = {"n": len(keep), "n_events": int(ev.sum()),
           "n_iatlas_samples_with_expression": int(n_ext),
           "event_source": "iAtlas OS_STATUS (true events; replaces the "
                           "0 < OS <= 400 d approximation)",
           "matching": "greedy one-to-one expression-profile correlation; "
                       f"min margin {min(r['margin'] for r in rows if r['margin'] is not None):.3f}",
           "os_agreement": f"{len(agree)}/{n_os} local OS_days agree with "
                           "iAtlas OS_MONTHS (±60 d)",
           "n_local_os_days_zero": int((a.obs['OS_days'].values <= 0).sum())}

    for name, score in (("impres", impres_score(Xk)), ("gep", gep_score(Xk))):
        z = (score - score.mean()) / score.std()
        df = pd.DataFrame({"dur": dur, "event": ev, "sig_z": z})
        cph = CoxPHFitter()
        cph.fit(df, duration_col="dur", event_col="event")
        s_ = cph.summary.loc["sig_z"]
        out[f"{name}_cox"] = {
            "HR_per_SD": round(float(np.exp(s_["coef"])), 3),
            "CI_low": round(float(np.exp(s_["coef lower 95%"])), 3),
            "CI_high": round(float(np.exp(s_["coef upper 95%"])), 3),
            "p": float(s_["p"])}
        hi = score > np.median(score)
        lr = logrank_test(dur[hi], dur[~hi], ev[hi], ev[~hi])
        out[f"{name}_km"] = {"p": float(lr.p_value)}
    resp = np.array([r["local_response"] for r in keep])
    lr = logrank_test(dur[resp == 1], dur[resp == 0], ev[resp == 1],
                      ev[resp == 0])
    out["recist_km"] = {"p": float(lr.p_value)}
    out["matches"] = rows

    with open(RES / "hugo_survival_true.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "matches"},
                     indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
