"""Stage-stratified Cox + landmark time-window analysis: the PH-violation
solution for the TCGA multivariate models (SKCM/BRCA stage covariate).

Analysis A (PH-robust primary sensitivity): stage-stratified Cox
  (signature z + age z; stage-specific baseline hazards, no PH assumption
  on stage), 12 signatures x 4 cancers, BH within cancer; compared against
  the archived standard multivariate results.
Analysis B (time-stability): 24-month landmark windows. Early window = all
  patients (events/censoring within 24 months). Late window = landmark at
  24 months (patients at risk at month 24, follow-up rebased). Standard
  multivariate (sig z + age z + stage z) per window; signature HRs per
  window show whether the signature effect drifts over time.

Output: data/tcga/stage_stratified_ph.json
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
TCGA = REPO / "data" / "tcga"
SPLIT_MONTHS = 24

spec = importlib.util.spec_from_file_location(
    "phc", REPO / "scripts" / "tcga_ph_check_schoenfeld.py")
phc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phc)  # loaders + SIGS + parse_stage
SIGS = phc.SIGS


def fit_std(d):
    cph = CoxPHFitter()
    cph.fit(d, duration_col="dur", event_col="event")
    s = cph.summary.loc["sig"]
    return {"HR": round(float(np.exp(s["coef"])), 3), "p": float(s["p"])}


def bh(ps):
    ps = np.asarray(ps, dtype=float)
    order = np.argsort(ps)
    m = len(ps)
    q = np.empty(m)
    q[order] = np.minimum.accumulate(
        (ps[order] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def main():
    out = {"split_months": SPLIT_MONTHS, "cancers": {}}
    for cancer in ["SKCM", "BRCA", "LUAD", "COAD"]:
        if cancer == "SKCM":
            df = phc.load_skcm()
        elif cancer == "BRCA":
            df = phc.load_brca_luad("BRCA")
        elif cancer == "LUAD":
            df = phc.load_brca_luad("LUAD")
        else:
            df = phc.load_coad()
        print(f"== {cancer}: n={len(df)} events={int(df['event'].sum())}",
              flush=True)
        # archived multivariate for comparison
        archived = json.load(open(TCGA / cancer / f"{cancer}_cox_multivariate.json",
                                  encoding="utf-8"))["results"]
        entry = {"n": int(len(df)), "events": int(df["event"].sum()),
                 "stratified": {}, "landmark": {}}
        # ---- Analysis A: stage-stratified ----
        hrs, pvals, names = [], [], []
        for sname, sgenes in SIGS.items():
            present = [g for g in sgenes if g in df.columns]
            if len(present) < 2:
                continue
            d = df[["dur", "event", "age", "stage"] + present].copy()
            d["sig"] = d[present].mean(axis=1)
            d = d.drop(columns=present)
            for c in ["sig", "age"]:
                d[c] = (d[c] - d[c].mean()) / d[c].std()
            try:
                cph = CoxPHFitter()
                cph.fit(d, duration_col="dur", event_col="event",
                        strata=["stage"])
                s = cph.summary.loc["sig"]
                hr = round(float(np.exp(s["coef"])), 3)
                pv = float(s["p"])
            except Exception as e:
                entry["stratified"][sname] = {"error": repr(e)[:80]}
                continue
            hrs.append(hr); pvals.append(pv); names.append(sname)
            entry["stratified"][sname] = {"HR": hr, "p": pv}
        qs = bh(pvals)
        n_sig = 0
        for name, hr, pv, q in zip(names, hrs, pvals, qs):
            entry["stratified"][name]["p_bh"] = float(q)
            if q < 0.05:
                n_sig += 1
            entry["stratified"][name]["archived_HR"] = archived.get(
                name, {}).get("HR")
            entry["stratified"][name]["archived_p_bh"] = archived.get(
                name, {}).get("p_bh")
        entry["stratified_n_sig"] = n_sig
        entry["archived_n_sig"] = sum(
            1 for v in archived.values()
            if isinstance(v, dict) and v.get("p_bh", 1) < 0.05)
        # ---- Analysis B: landmark windows ----
        early = df[df["dur"] <= SPLIT_MONTHS]
        late_all = df[df["dur"] > SPLIT_MONTHS].copy()
        late_all["dur"] = late_all["dur"] - SPLIT_MONTHS
        for sname, sgenes in SIGS.items():
            present = [g for g in sgenes if g in df.columns]
            if len(present) < 2:
                continue
            windows = {}
            for label, sub in (("early", early), ("late", late_all)):
                d = sub[["dur", "event", "age", "stage"] + present].copy()
                d["sig"] = d[present].mean(axis=1)
                d = d.drop(columns=present)
                if d["event"].sum() < 5 or len(d) < 30:
                    windows[label] = {"skipped": "too few events/n"}
                    continue
                for c in ["sig", "age", "stage"]:
                    d[c] = (d[c] - d[c].mean()) / d[c].std()
                try:
                    windows[label] = fit_std(d)
                    windows[label]["n"] = int(len(d))
                    windows[label]["events"] = int(d["event"].sum())
                except Exception as e:
                    windows[label] = {"error": repr(e)[:80]}
            if "early" in windows and "late" in windows and \
                    "HR" in windows.get("early", {}) and \
                    "HR" in windows.get("late", {}):
                windows["direction_stable"] = (
                    (windows["early"]["HR"] < 1) == (windows["late"]["HR"] < 1))
            entry["landmark"][sname] = windows
        out["cancers"][cancer] = entry
        print(f"   stratified BH-sig: {n_sig} (archived "
              f"{entry['archived_n_sig']})", flush=True)
    dst = TCGA / "stage_stratified_ph.json"
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written:", dst, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
