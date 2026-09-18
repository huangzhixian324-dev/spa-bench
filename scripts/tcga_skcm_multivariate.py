"""Dry experiment C: TCGA-SKCM multivariate Cox sensitivity (signature +
age + AJCC stage), addressing the "univariate confounding" attack on the
immune-signature prognostic gradient.

Data: refetched expression (44 signature genes, API) + refetch OS + local
pan-cancer clinical (AGE, AJCC_PATHOLOGIC_TUMOR_STAGE, 441/441 patients).
Stage parsed to ordinal 1-4 (leading roman numeral; unknown dropped).
Signatures and covariates z-standardized. BH within the 12 signatures.

Output: data/tcga/SKCM/SKCM_cox_multivariate.json
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

REPO = Path(__file__).resolve().parents[1]
TCGA = REPO / "data" / "tcga"
REFETCH = TCGA / "SKCM" / "refetch"

SIGS = {
    'IMPRES-like': ['PDCD1','CD27','CTLA4','CD28','CD86','CD80','CD274','LAG3','HAVCR2','TIGIT','TNFRSF9','ICOS'],
    'CD8 T cells': ['CD8A','CD8B'], 'B cells': ['CD19','CD79A','MS4A1'],
    'NK cells': ['NKG7','KLRD1','KLRF1'], 'IFNG response': ['IFNG','CXCL10','CXCL9','IDO1','STAT1'],
    'Cytolytic': ['GZMA','PRF1','GNLY'], 'Checkpoint': ['PDCD1','CTLA4','LAG3','TIGIT','HAVCR2'],
    'TLS': ['CXCL13','CCL19','CCL21'], 'Treg': ['FOXP3','IL2RA'],
    'M1 Macrophage': ['NOS2','IL12A','TNF'], 'Myeloid': ['CD14','CD33','ITGAM'],
    'Stromal': ['COL1A1','COL1A2','COL3A1','ACTA2','FAP'],
}
STAGE_MAP = {"I": 1, "II": 2, "III": 3, "IV": 4}


def parse_stage(v):
    if not isinstance(v, str):
        return np.nan
    m = re.match(r"(?i)stage\s+([IV]+)", v.strip())
    if not m:
        return np.nan
    return STAGE_MAP.get(m.group(1).upper(), np.nan)


def main():
    expr = pd.read_csv(REFETCH / "expression.tsv.gz", sep="\t", index_col=0)
    os_cli = pd.read_csv(REFETCH / "clinical.tsv", sep="\t", index_col=0)
    local = pd.read_csv(TCGA / "clinical_skcm.txt", sep="\t", comment="#",
                        index_col=0)

    # sample -> patient; first tumour sample per patient
    pats = ["-".join(c.split("-")[:3]) for c in expr.columns]
    seen = {}
    for c, p in zip(expr.columns, pats):
        if p not in seen:
            seen[p] = c
    expr2 = expr[list(seen.values())]
    expr2.columns = list(seen.keys())
    df = expr2.T.join(os_cli[["OS_STATUS", "OS_MONTHS"]], how="inner")
    df = df.join(local[["AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE"]], how="left")
    df["event"] = df["OS_STATUS"].astype(str).str.upper() \
        .str.contains("DECEASED").astype(int)
    df["dur"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
    df["age"] = pd.to_numeric(df["AGE"], errors="coerce")
    df["stage"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(parse_stage)
    df = df.dropna(subset=["dur", "age", "stage"])
    df = df[df["dur"] > 0]
    print(f"n={len(df)} events={int(df['event'].sum())} "
          f"stage dist={df['stage'].value_counts().sort_index().to_dict()}",
          flush=True)

    genes = [g for g in expr.index if g in df.columns] if False else None
    results = {}
    for sname, sgenes in SIGS.items():
        present = [g for g in sgenes if g in df.columns]
        if len(present) < 2:
            continue
        d = df[["dur", "event", "age", "stage"] + present].copy()
        d["sig"] = d[present].mean(axis=1)
        d = d.drop(columns=present)
        # z-standardize
        for c in ["sig", "age", "stage"]:
            d[c] = (d[c] - d[c].mean()) / d[c].std()
        cph = CoxPHFitter()
        cph.fit(d, duration_col="dur", event_col="event")
        s = cph.summary.loc["sig"]
        results[sname] = {
            "HR": round(float(np.exp(s["coef"])), 3),
            "p": float(s["p"]),
            "HR_age": round(float(np.exp(cph.params_["age"])), 3),
            "p_age": float(cph.summary.loc["age", "p"]),
            "HR_stage": round(float(np.exp(cph.params_["stage"])), 3),
            "p_stage": float(cph.summary.loc["stage", "p"]),
        }
    tested = list(results)
    ps = np.array([results[s]["p"] for s in tested])
    order = np.argsort(ps)
    m = len(ps)
    bh = np.empty(m)
    bh[order] = np.minimum.accumulate(
        (ps[order] * m / (np.arange(m) + 1))[::-1])[::-1]
    for s, q in zip(tested, bh):
        results[s]["p_bh"] = float(min(q, 1.0))
    n_sig = sum(1 for v in results.values() if v["p_bh"] < 0.05)
    hrs = sorted(v["HR"] for v in results.values() if v["p_bh"] < 0.05)
    print(f"multivariate (sig + age + stage): {n_sig}/{len(results)} "
          f"BH-significant; HR range {hrs[0]:.3f}-{hrs[-1]:.3f}"
          if hrs else "no BH-significant signature", flush=True)

    out = {"cohort": "TCGA-SKCM",
           "model": "multivariate Cox: signature (z) + AGE (z) + AJCC stage "
                    "(ordinal 1-4, z); BH within 12 signatures",
           "n": int(len(df)), "events": int(df["event"].sum()),
           "note": "addresses the univariate-confounding attack on the "
                   "gradient; complements the sample-type decomposition",
           "results": results}
    dst = TCGA / "SKCM" / "SKCM_cox_multivariate.json"
    json.dump(out, open(dst, "w"), indent=1)
    print("WROTE", dst, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
