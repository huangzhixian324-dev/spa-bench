"""Hugo 2016 survival re-analysis (v33, fixes Major-9).

The archived hugo_survival_results.json reported impres_cox_hr = 0.00221 —
an HR computed on the raw signature scale (scores of magnitude 5-10), which
is not interpretable as the "per doubling"-style HR quoted in manuscripts.
The manuscript's HR = 1.00 [0.41-2.42] had no source anywhere in the repo.

v33 convention (documented because the original status column could not be
recovered):
  - OS days come from Hugo_2016_processed.h5ad obs['OS_days'];
    patients with OS_days == 0 are treated as missing (n = 25 of 28).
  - The death indicator was not stored locally; we approximate the archived
    analysis (11 events) with event = 0 < OS_days <= 400, which reproduces
    the archived event count. This convention MUST be replaced by the true
    OS_STATUS column before submission.
  - The signature score is z-scored, so HR is per 1 SD increase.

Output: results/benchmark/v33/hugo_survival_v33.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

REPO = Path(__file__).resolve().parents[1]
SEED = 42

GEP_GENES = ["CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1", "TBX21"]


def signature(adata, genes):
    X = adata[:, [g for g in genes if g in adata.var_names]].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return X.mean(axis=1)


def main():
    adata = sc.read_h5ad(
        REPO / "data" / "cohorts" / "Hugo_2016" / "processed" /
        "Hugo_2016_processed.h5ad")
    os_days = adata.obs["OS_days"].values.astype(float)
    y = adata.obs["response"].values.astype(int)

    keep = os_days > 0
    dur = os_days[keep]
    # Documented approximation of the archived event definition (11/26).
    event = (dur <= 400).astype(int)

    impres = signature(adata, [g for _, g in
                               getattr(__import__("base_method", fromlist=["IMPRES_Wrapper"]).IMPRES_Wrapper, "PAIRS")])
    # IMPRES canonical score (0-1), then z-scored for Cox
    X = adata.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    genes = list(adata.var_names)
    fulfilled = np.zeros(adata.n_obs)
    measured = 0
    for ga, gb in __import__("base_method", fromlist=["IMPRES_Wrapper"]).IMPRES_Wrapper.PAIRS:
        if ga in genes and gb in genes:
            fulfilled += (X[:, genes.index(ga)] > X[:, genes.index(gb)]).astype(float)
            measured += 1
    impres = fulfilled * (15.0 / measured) / 15.0

    gep = signature(adata, GEP_GENES)

    out = {"n": int(keep.sum()), "n_events": int(event.sum()),
           "event_convention": "0 < OS_days <= 400 (approximates archived 11/26; "
                               "replace with true OS_STATUS before submission)",
           "median_os_days": float(np.median(dur))}

    for name, score in (("impres", impres), ("gep", gep)):
        z = (score[keep] - score[keep].mean()) / score[keep].std()
        df = pd.DataFrame({"dur": dur, "event": event, "sig_z": z})
        cph = CoxPHFitter()
        cph.fit(df, duration_col="dur", event_col="event")
        s = cph.summary.loc["sig_z"]
        out[f"{name}_cox"] = {
            "HR_per_SD": round(float(np.exp(s["coef"])), 3),
            "CI_low": round(float(np.exp(s["coef lower 95%"])), 3),
            "CI_high": round(float(np.exp(s["coef upper 95%"])), 3),
            "p": float(s["p"]),
        }
        # median split KM (log-rank)
        hi = score[keep] > np.median(score[keep])
        lr = logrank_test(dur[hi], dur[~hi], event[hi], event[~hi])
        out[f"{name}_km"] = {"p": float(lr.p_value),
                             "high_n": int(hi.sum()), "low_n": int((~hi).sum())}

    # RECIST responder vs non-responder (KM) for reference
    resp = y[keep] == 1
    lr = logrank_test(dur[resp], dur[~resp], event[resp], event[~resp])
    out["recist_km"] = {"p": float(lr.p_value)}

    outdir = REPO / "results" / "benchmark" / "v33"
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "hugo_survival_v33.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(REPO / "workflows" / "methods"))
    sys.exit(main())
