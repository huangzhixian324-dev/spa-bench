"""Gide 2019 clinical-utility metrics for the FDR-significant cells
(reviewer item R2-M1): OS stratification (median-split log-rank + Cox HR
per SD) and decision-analytic metrics (Brier, DCA net benefit, IDI) using
OOF predictions regenerated with the pipeline's own functions
(rerun_v33.load_cohort / cv_predict / nested_cv_predict, frozen modal
hyperparameters from tune_primary).

OS source: cBioPortal API mel_iatlas_gide_2019 patient clinical data
(OS_STATUS, OS_MONTHS; fetched 2026-09-27). Survival metrics are computed
only for patients matched to the 73 archived samples.

Output: results/benchmark/v33/gide_utility_metrics.json
"""
import importlib.util
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from sklearn.metrics import brier_score_loss, roc_auc_score

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
spec = importlib.util.spec_from_file_location(
    "rerun_v33", REPO / "scripts" / "rerun_v33.py")
rerun = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rerun)
import scanpy as sc  # noqa: E402  (after rerun import sets sys.path)

OUT = rerun.OUT / "gide_utility_metrics.json"
BASE = "https://www.cbioportal.org/api"
THRESHOLDS = [0.10, 0.20, 0.30, 0.50]


def fetch_os():
    urllib.request.install_opener(urllib.request.build_opener(
        urllib.request.ProxyHandler({})))
    rows, page = {}, 0
    while True:
        url = (f"{BASE}/studies/mel_iatlas_gide_2019/clinical-data"
               f"?clinicalDataType=PATIENT&projection=DETAILED"
               f"&pageSize=100000&page={page}")
        with urllib.request.urlopen(url, timeout=60) as r:
            d = json.load(r)
        if not d:
            break
        for x in d:
            pid = x["patientId"]
            rows.setdefault(pid, {})[x["clinicalAttributeId"]] = x.get("value")
        if len(d) < 100000:
            break
        page += 1
    os_rows = {}
    for pid, attrs in rows.items():
        st = attrs.get("OS_STATUS")
        mo = attrs.get("OS_MONTHS")
        if st is not None and mo is not None:
            os_rows[pid] = {"OS_STATUS": st, "OS_MONTHS": mo}
    return os_rows


def main():
    adata, folds, y = rerun.load_cohort("Gide_2019_cBio")
    adata.obs["response"] = y
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    var_names = list(adata.var_names)
    sample_ids = list(adata.obs_names)
    patients = [s.split("_")[0] for s in sample_ids]

    frozen = {}
    for m in ("ElasticNet", "ElasticNet_Var"):
        frozen[m] = rerun.tune_primary(adata, folds, m)
    print("frozen params:", frozen, flush=True)

    scores = {}
    for m in ("GEP", "TIDE", "IMPRES"):
        yt, yp = rerun.cv_predict(adata, folds, m)
        scores[m] = (yt, yp)
    for m in ("ElasticNet", "ElasticNet_Var"):
        yt, yp, _ = rerun.nested_cv_predict(adata, folds, m,
                                            params=frozen[m])
        scores[m] = (yt, yp)
    # PD-L1 proxy: raw CD274 expression (higher = more positive)
    if "CD274" in var_names:
        cd = X[:, var_names.index("CD274")]
        scores["PD-L1 (CD274)"] = (y.astype(int), cd)

    os_rows = fetch_os()
    print("OS matched patients:", len(os_rows), flush=True)

    os_map, keep_idx = {}, []
    for i, pat in enumerate(patients):
        if pat in os_rows:
            keep_idx.append(i)
            st = os_rows[pat]["OS_STATUS"]
            os_map[i] = (1 if "DECEASED" in str(st).upper() else 0,
                         float(os_rows[pat]["OS_MONTHS"]))
    df_os = pd.DataFrame({
        "patient": [patients[i] for i in keep_idx],
        "event": [os_map[i][0] for i in keep_idx],
        "os_months": [os_map[i][1] for i in keep_idx],
    }).set_index("patient")
    print(f"OS-matched samples: {len(keep_idx)} "
          f"events={int(df_os['event'].sum())}", flush=True)

    z = lambda v: (np.asarray(v, dtype=float) - np.mean(v)) / np.std(v)

    methods, out_metrics = {}, {}
    for m, (yt, yp) in scores.items():
        yp = np.asarray(yp, dtype=float)
        yv = np.asarray(yt, dtype=int)
        rec = {"auroc": round(float(roc_auc_score(yv, yp)), 4)}
        if len(keep_idx) == len(yp):
            # survival stratification: median-split score
            med = np.median(yp)
            hi = yp >= med
            e = df_os["event"].values
            d = df_os["os_months"].values
            lr = logrank_test(d[hi], d[~hi], event_observed_A=e[hi],
                              event_observed_B=e[~hi])
            dd = pd.DataFrame({"dur": d, "event": e, "score": z(yp)})
            cph = CoxPHFitter()
            cph.fit(dd, duration_col="dur", event_col="event")
            s = cph.summary.loc["score"]
            rec["os_logrank_p"] = float(lr.p_value)
            rec["os_cox_HR_per_SD"] = round(float(np.exp(s["coef"])), 3)
            rec["os_cox_p"] = float(s["p"])
            rec["os_split"] = "median score"
            # decision-analytic metrics (probabilities only for ElasticNets)
            if m in ("ElasticNet", "ElasticNet_Var"):
                rec["brier"] = round(float(brier_score_loss(yv, yp)), 4)
                nb = {}
                for pt in THRESHOLDS:
                    pred = (yp >= pt).astype(int)
                    tp = int(((pred == 1) & (yv == 1)).sum())
                    fp = int(((pred == 1) & (yv == 0)).sum())
                    nb[str(pt)] = round(tp / len(yv) -
                                        fp / len(yv) * (pt / (1 - pt)), 4)
                treat_all = round(float(yv.mean()), 4)
                rec["dca_net_benefit"] = nb
                rec["dca_treat_all"] = treat_all
        out_metrics[m] = rec
        methods[m] = {"frozen_params": frozen.get(m)}
    # IDI (mean prob difference) for IMPRES/EN pairs with [0,1] scores
    idi = {}
    for m in ("IMPRES", "ElasticNet", "ElasticNet_Var"):
        yt, yp = scores[m]
        yp = np.asarray(yp, dtype=float)
        if yp.min() < 0 or yp.max() > 1:
            continue
        yv = np.asarray(yt, dtype=int)
        idi[m] = {"IDI": round(float(yp[yv == 1].mean() - yp[yv == 0].mean()), 4)}
    if "IMPRES" in idi:
        for m in ("ElasticNet", "ElasticNet_Var"):
            if m in idi and "IMPRES" in idi:
                idi[f"{m}_vs_IMPRES_NRI-proxy"] = round(
                    idi[m]["IDI"] - idi["IMPRES"]["IDI"], 4)

    out = {"cohort": "Gide_2019_cBio", "endpoint": "RECIST",
           "os_source": ("cBioPortal API mel_iatlas_gide_2019 patient "
                         "clinical data, fetched 2026-09-27"),
           "n_samples": int(len(y)), "n_os_matched": int(len(keep_idx)),
           "os_events": int(df_os["event"].sum()),
           "frozen_params": frozen,
           "metrics": out_metrics, "IDI": idi,
           "notes": ("Brier/DCA computed only for ElasticNet variants "
                     "(calibrated probabilities); fixed-scorer scores are "
                     "uncalibrated so only survival metrics are reported. "
                     "OS stratification uses a median score split."),
           }
    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written:", OUT, flush=True)
    print(json.dumps(out_metrics, indent=1)[:1800], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
