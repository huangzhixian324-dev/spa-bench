"""Generalized TCGA multivariate Cox (signature + age + stage) for BRCA/LUAD/COAD.
Mirrors scripts/tcga_skcm_multivariate.py protocol exactly:
- cBioPortal raw expression (BRCA/LUAD) or GDC STAR (COAD, if present)
- patient-level first tumor sample; OS from patient clinical; AGE + AJCC stage
- signature (z) + AGE (z) + stage (ordinal 1-4, z); BH within 12 signatures
Outputs: data/tcga/<C>/<C>_cox_multivariate.json
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
TCGA = REPO / "data" / "tcga"

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


def load_expr(clinical_df):
    """Return expr DataFrame (genes x patients) from the available source."""
    gdc = TCGA / "COAD" / "raw" / "data_mrna_gdc_star_fpkm.txt"
    return None, None  # placeholder (unused for BRCA/LUAD)


def run(cancer, expr_path, sep="\t", is_gdc=False):
    raw = TCGA / cancer / "raw"
    expr = pd.read_csv(expr_path, sep=sep, index_col=0)
    # drop cBioPortal hugo placeholder rows if present
    expr = expr[~expr.index.astype(str).str.contains("ESR1|?|", regex=False)]
    if not is_gdc:
        expr = expr[~expr.index.astype(str).str.startswith("?")]
    patients = ["-".join(c.split("-")[:3]) for c in expr.columns]
    seen = {}
    for col, p in zip(expr.columns, patients):
        if p not in seen:
            seen[p] = col
    expr2 = expr[list(seen.values())]
    expr2.columns = list(seen.keys())

    cli = pd.read_csv(raw / "data_clinical_patient.txt", sep="\t",
                      comment="#", index_col=0, dtype=str)
    df = expr2.T.join(cli[["OS_STATUS", "OS_MONTHS", "AGE",
                           "AJCC_PATHOLOGIC_TUMOR_STAGE"]], how="inner")
    df["event"] = df["OS_STATUS"].astype(str).str.upper() \
        .str.contains("DECEASED").astype(int)
    df["dur"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
    df["age"] = pd.to_numeric(df["AGE"], errors="coerce")
    df["stage"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(parse_stage)
    df = df.dropna(subset=["dur", "age", "stage"])
    df = df[df["dur"] > 0]
    print(f"{cancer}: n={len(df)} events={int(df['event'].sum())}", flush=True)

    results = {}
    for sname, sgenes in SIGS.items():
        present = [g for g in sgenes if g in df.columns]
        if len(present) < 2:
            continue
        d = df[["dur", "event", "age", "stage"] + present].copy()
        d["sig"] = d[present].mean(axis=1)
        d = d.drop(columns=present)
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
    rng = f"{hrs[0]:.3f}-{hrs[-1]:.3f}" if hrs else "none"
    print(f"{cancer} multivariate: {n_sig}/{len(results)} BH-sig; HR {rng}",
          flush=True)
    out = {"cohort": f"TCGA-{cancer}",
           "model": "multivariate Cox: signature (z) + AGE (z) + AJCC stage "
                    "(ordinal 1-4, z); BH within 12 signatures",
           "n": len(df), "events": int(df["event"].sum()),
           "note": ("extends the SKCM multivariate sensitivity to this cancer "
                    "type; addresses the univariate-confounding attack on the "
                    "gradient"),
           "results": results}
    (TCGA / cancer / f"{cancer}_cox_multivariate.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    return out


def main():
    summary = {}
    for cancer in ["BRCA", "LUAD"]:
        expr_path = TCGA / cancer / "raw" / "data_mrna_seq_v2_rsem.txt"
        summary[cancer] = run(cancer, expr_path)
    coad_gdc = TCGA / "COAD" / "raw" / "data_mrna_gdc_star_fpkm.txt"
    if coad_gdc.exists():
        summary["COAD"] = run("COAD", coad_gdc, is_gdc=True)
    else:
        summary["COAD"] = {"status": "pending GDC refetch "
                                     "(data_mrna_gdc_star_fpkm.txt absent)"}
        print("COAD: GDC expression absent — multivariate pending refetch",
              flush=True)
    (TCGA / "cox_multivariate_all.json").write_text(
        json.dumps(summary, indent=1), encoding="utf-8")
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
