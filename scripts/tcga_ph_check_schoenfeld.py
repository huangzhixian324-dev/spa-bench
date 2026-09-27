"""Schoenfeld-residual proportional-hazards check for the TCGA multivariate
Cox models reported in Supplementary Table S15 (SKCM/BRCA/LUAD/COAD).

Rebuilds each cancer's model frame with the EXACT construction of the
archived scripts (tcga_skcm_multivariate.py / tcga_multivariate_all.py /
tcga_coad_multivariate_gdc.py), refits each of the 12 signature models
identically (CoxPHFitter, z-standardized sig/age/stage), then runs
lifelines.statistics.proportional_hazard_test (rank transform) per model.

Output: data/tcga/ph_schoenfeld_check.json
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
TCGA = REPO / "data" / "tcga"
RAW = TCGA / "COAD" / "raw"

spec = importlib.util.spec_from_file_location(
    "tma", REPO / "scripts" / "tcga_multivariate_all.py")
tma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tma)  # SIGS + parse_stage (identical in all scripts)
SIGS, parse_stage = tma.SIGS, tma.parse_stage


def frame_from_expr(expr, cli, cols):
    patients = ["-".join(c.split("-")[:3]) for c in expr.columns]
    seen = {}
    for col, p in zip(expr.columns, patients):
        if p not in seen:
            seen[p] = col
    expr2 = expr[list(seen.values())]
    expr2.columns = list(seen.keys())
    df = expr2.T.join(cli[cols], how="inner")
    df["event"] = df["OS_STATUS"].astype(str).str.upper() \
        .str.contains("DECEASED").astype(int)
    df["dur"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
    df["age"] = pd.to_numeric(df["AGE"], errors="coerce")
    df["stage"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(parse_stage)
    df = df.dropna(subset=["dur", "age", "stage"])
    return df[df["dur"] > 0]


def load_brca_luad(cancer):
    raw = TCGA / cancer / "raw"
    expr = pd.read_csv(raw / "data_mrna_seq_v2_rsem.txt", sep="\t", index_col=0)
    expr = expr[~expr.index.astype(str).str.contains("ESR1|?|", regex=False)]
    expr = expr[~expr.index.astype(str).str.startswith("?")]
    cli = pd.read_csv(raw / "data_clinical_patient.txt", sep="\t",
                      comment="#", index_col=0, dtype=str)
    return frame_from_expr(expr, cli, ["OS_STATUS", "OS_MONTHS", "AGE",
                                       "AJCC_PATHOLOGIC_TUMOR_STAGE"])


def load_skcm():
    refetch = TCGA / "SKCM" / "refetch"
    expr = pd.read_csv(refetch / "expression.tsv.gz", sep="\t", index_col=0)
    os_cli = pd.read_csv(refetch / "clinical.tsv", sep="\t", index_col=0)
    local = pd.read_csv(TCGA / "clinical_skcm.txt", sep="\t", comment="#",
                        index_col=0)
    patients = ["-".join(c.split("-")[:3]) for c in expr.columns]
    seen = {}
    for c, p in zip(expr.columns, patients):
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
    return df[df["dur"] > 0]


def load_coad():
    expr = pd.read_csv(RAW / "data_mrna_gdc_star_fpkm.txt", sep="\t",
                       index_col=0)
    expr = expr[~expr.index.astype(str).str.startswith("?")]
    patients = ["-".join(c.split("-")[:3]) for c in expr.columns]
    seen = {}
    for col, p in zip(expr.columns, patients):
        if p not in seen:
            seen[p] = col
    expr2 = expr[list(seen.values())]
    expr2.columns = list(seen.keys())
    side = pd.read_csv(RAW / "data_clinical_patient_gdc_multiv.txt", sep="\t",
                       dtype=str, index_col=0)
    cli = side[["OS_STATUS", "OS_MONTHS", "AGE",
                "AJCC_PATHOLOGIC_TUMOR_STAGE"]]
    df = expr2.T.join(cli, how="inner")
    df["event"] = df["OS_STATUS"].astype(str).str.upper() \
        .str.contains("DECEASED").astype(int)
    df["dur"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
    df["age"] = pd.to_numeric(df["AGE"], errors="coerce")
    df["stage"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(parse_stage)
    df = df.dropna(subset=["dur", "age", "stage"])
    return df[df["dur"] > 0]


LOADERS = {"SKCM": load_skcm, "BRCA": load_brca_luad, "LUAD": load_brca_luad,
           "COAD": load_coad}


def main():
    out = {"model": "multivariate Cox: signature (z) + AGE (z) + AJCC stage "
                    "(z); PH test = lifelines proportional_hazard_test, "
                    "rank transform, per covariate (sig/age/stage)",
           "alpha": 0.05, "cancers": {}}
    for cancer in ["SKCM", "BRCA", "LUAD", "COAD"]:
        if cancer == "BRCA":
            df = load_brca_luad("BRCA")
        elif cancer == "LUAD":
            df = load_brca_luad("LUAD")
        else:
            df = LOADERS[cancer]()
        print(f"== {cancer}: n={len(df)} events={int(df['event'].sum())}",
              flush=True)
        entry = {"n": int(len(df)), "events": int(df["event"].sum()),
                 "signatures": {}}
        n_violation_models = 0
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
            ph = proportional_hazard_test(cph, d, time_transform="rank")
            per_cov = {k: float(v) for k, v in ph.summary["p"].items()}
            violated = [k for k, v in per_cov.items() if v < 0.05]
            if violated:
                n_violation_models += 1
            entry["signatures"][sname] = {
                "ph_p": {k: round(v, 4) for k, v in per_cov.items()},
                "violated_at_0.05": violated,
            }
        entry["models_with_any_violation"] = n_violation_models
        entry["models_total"] = len(entry["signatures"])
        out["cancers"][cancer] = entry
        print(f"   PH violations (p<0.05): {n_violation_models}/"
              f"{len(entry['signatures'])} models", flush=True)
    dst = TCGA / "ph_schoenfeld_check.json"
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written:", dst, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
