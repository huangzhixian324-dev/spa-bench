"""CDS v2.0 hypothesis-free validation matrix — two phases
(reviewer items R1-M4/R2-M4/R3-M2, strengthening round).

Phase 1 (circular-by-construction, 4 endpoints): cytolytic surrogate
  = median split of GZMA/PRF1 mean, on Liu 2019 (mel_iatlas_liu_2019,
  cBioPortal), IMvigor210 (blca_iatlas_imvigor210_2017, cBioPortal),
  STAD_PRJEB25780 (local HGNC h5ad), GSE274975 (local h5ad).
Phase 2 (clinical RECIST controls, 4 endpoints): the same cohorts'
  RECIST-style responder labels.

Per endpoint: CDS v2.0.0 hypothesis-free statistic (p90 of genome-wide MI,
200 label permutations, top-2,000-variance features) + v1 composite +
declared-gene specificity. 200 permutations per endpoint.

Output: results/benchmark/v33/cds_v2_hypfree_matrix.json (checkpointed)
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
sys.path.insert(0, str(REPO / "cds_tool"))
import cds as cds_mod  # noqa: E402

PERM = 200
EP = ["GZMA", "PRF1"]
DROP_VALUES = {"", "NA", "UNKNOWN", "NAN", "NONE"}
POS_VALUES = {"1", "PR", "CR", "RESPONDER", "TRUE", "R", "PARTIAL RESPONSE", "COMPLETE RESPONSE"}
OUT = REPO / "results" / "benchmark" / "v33" / "cds_v2_hypfree_matrix.json"
CK = OUT.with_suffix(".checkpoint.json")


def save(out):
    CK.write_text(json.dumps(out, indent=1), encoding="utf-8")
    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")

# ---- cBioPortal fetchers reused from the E3 literature backtest ----
spec = importlib.util.spec_from_file_location(
    "e3", REPO / "scripts" / "e3_literature_backtest.py")
e3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e3)

STAD_P = REPO / "data/cohorts/STAD_PRJEB25780/processed/STAD_PRJEB25780_HGNC.h5ad"
GSE_P = REPO / "data/cohorts/GSE274975/processed/GSE274975_HGNC.h5ad"
RESP_CANDIDATES = ["RECIST", "RESPONDER", "response", "CLINICAL_BENEFIT",
                   "RECIST_RESPONDER", "BOR"]

RECIST_ATTR_FOUND = {}


def pick_responder(attrs):
    for cand in RESP_CANDIDATES:
        for k, v in attrs.items():
            if cand.lower() in k.lower():
                RECIST_ATTR_FOUND[cand] = k
                return v
    return None


def read_expr(study):
    cache = REPO / "data" / "external" / study / "expression.tsv.gz"
    df = pd.read_csv(cache, sep="	", index_col=0)
    if df.shape[1] == 0:
        df = pd.read_csv(cache, index_col=0)
    df = df.dropna(axis=0)
    return df


def cyt_label(X, genes):
    gidx = [genes.index(g) for g in EP if g in genes]
    cyt = X[:, gidx].mean(axis=1)
    return (cyt > np.median(cyt)).astype(int)


def run_endpoint(name, X, y, genes, out):
    genes_l = [str(g) for g in genes]
    print(f"== {name}: {X.shape[0]} x {X.shape[1]} pos={int(np.sum(y))}",
          flush=True)
    r = cds_mod.circularity_detection_score(
        X, np.asarray(y).astype(int), EP, genes_l, permutation=PERM,
        hypothesis_free=True)
    v1 = {"CDS": r["CDS"], "risk": r["risk_level"]}
    v2 = r.get("v2", {})
    rec = {"v1": v1,
           "hyp_free_p": v2.get("hypothesis_free", {}).get("empirical_p"),
           "declared_pct": v2.get("declared_gene_specificity", {}).get(
               "mean_percentile"),
           "declared_p": v2.get("declared_gene_specificity", {}).get(
               "empirical_p"),
           "perm_c3_p": v2.get("permutation_c3", {}).get("empirical_p")}
    out["endpoints"][name] = rec
    print(f"   v1 CDS={rec['v1']['CDS']} [{rec['v1']['risk']}] | "
          f"hyp-free p={rec['hyp_free_p']} | declared pct={rec['declared_pct']} "
          f"(p={rec['declared_p']})", flush=True)
    return rec


def main():
    out = {"meta": {"cds_version": cds_mod.__version__,
                    "permutation": PERM,
                    "statistic": "hypothesis-free p90 genome-wide MI, "
                                 "top-2,000-variance features",
                    "phases": {"1": "circular cytolytic surrogates (4)",
                               "2": "clinical RECIST controls (4)"}},
           "endpoints": {}}
    if CK.exists():
        out = json.load(open(CK, encoding="utf-8"))
        print("resumed from checkpoint", flush=True)

    # ---------- Phase 1: circular ----------
    if "Liu_cytolytic" not in out["endpoints"]:
        expr = read_expr("mel_iatlas_liu_2019")
        clin = e3.fetch_clinical("mel_iatlas_liu_2019")
        genes = [str(g) for g in expr.index]
        X = expr.T.values
        y = cyt_label(X, genes)
        run_endpoint("Liu_cytolytic", X, y, genes, out)
        CK.write_text(json.dumps(out, indent=1), encoding="utf-8")
    if "IMvigor_cytolytic" not in out["endpoints"]:
        expr = read_expr("blca_iatlas_imvigor210_2017")
        clin = e3.fetch_clinical("blca_iatlas_imvigor210_2017")
        genes = [str(g) for g in expr.index]
        X = expr.T.values
        y = cyt_label(X, genes)
        run_endpoint("IMvigor_cytolytic", X, y, genes, out)
        CK.write_text(json.dumps(out, indent=1), encoding="utf-8")
    if "STAD_cytolytic" not in out["endpoints"]:
        a = sc.read_h5ad(STAD_P)
        genes = [str(g) for g in a.var_names]
        X = (a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X))
        y = cyt_label(X, genes)
        run_endpoint("STAD_cytolytic", X, y, genes, out)
        CK.write_text(json.dumps(out, indent=1), encoding="utf-8")
    if "GSE274975_cytolytic" not in out["endpoints"]:
        a = sc.read_h5ad(GSE_P)
        genes = [str(g) for g in a.var_names]
        X = (a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X))
        y = cyt_label(X, genes)
        run_endpoint("GSE274975_cytolytic", X, y, genes, out)
        CK.write_text(json.dumps(out, indent=1), encoding="utf-8")

    # ---------- Phase 2: clinical RECIST controls ----------
    for study, name in (("mel_iatlas_liu_2019", "Liu_RECIST"),
                        ("blca_iatlas_imvigor210_2017", "IMvigor_RECIST")):
        if name in out["endpoints"]:
            continue
        # responder labels: local iAtlas clinical export first (RECIST
        # categories live there for Liu; the API omits them), API fallback
        expr = read_expr(study)
        genes = [str(g) for g in expr.index]
        resp_series, best = None, None
        local_ts = REPO / "data" / "external" / study / "clinical.tsv"
        if local_ts.exists():
            cl = pd.read_csv(local_ts, sep="	")
            for cand in ("RESPONSE", "RECIST", "RESPONDER",
                         "CLINICAL_BENEFIT"):
                if cand in cl.columns:
                    resp_series = cl.set_index("SAMPLE_ID")[cand]
                    best = cand
                    break
        if resp_series is None:
            clin = e3.fetch_clinical(study)
            best, hits = pick_responder(clin)
            print(f"  {name}: API responder attribute {best} ({hits})",
                  flush=True)
            if best:
                resp_series = pd.Series(
                    {sid: attrs.get(best) for sid, attrs in clin.items()})
        if resp_series is None:
            print(f"  {name}: no responder attribute found — skipped",
                  flush=True)
            continue
        print(f"  {name}: responder column = {best}", flush=True)
        y, keep, vals = [], [], {}
        for c in map(str, expr.columns):
            v = resp_series.get(c)
            v = "" if v is None else str(v).strip().upper()
            vals[v] = vals.get(v, 0) + 1
            if v in DROP_VALUES or v == "NAN":
                continue
            keep.append(c)
            y.append(1 if v in POS_VALUES else 0)
        print(f"  {name}: value distribution {vals}; kept {len(keep)}",
              flush=True)
        if len(keep) < 30:
            print(f"  {name}: too few annotated samples — skipped",
                  flush=True)
            continue
        Xsub = expr[keep].T.values
        run_endpoint(name, Xsub, np.array(y), genes, out)
        out["endpoints"][name]["responder_column"] = best
        out["endpoints"][name]["value_distribution"] = vals
        save(out)
    for local, name in ((STAD_P, "STAD_RECIST"), (GSE_P, "GSE274975_RECIST")):
        if name in out["endpoints"]:
            continue
        a = sc.read_h5ad(local)
        genes = [str(g) for g in a.var_names]
        X = (a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X))
        y = a.obs["response"].values.astype(int)
        run_endpoint(name, X, y, genes, out)
        CK.write_text(json.dumps(out, indent=1), encoding="utf-8")

    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written:", OUT, flush=True)
    print("RECIST attribute map:", RECIST_ATTR_FOUND, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
