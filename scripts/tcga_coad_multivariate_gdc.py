"""COAD multivariate Cox (signature + age + stage) from the refetched GDC matrix.

Companion to scripts/tcga_multivariate_all.py for the GDC STAR source. The
cBioPortal-format clinical file (AGE / AJCC_PATHOLOGIC_TUMOR_STAGE columns)
does not exist for COAD, so age and stage are fetched from the GDC cases
endpoint (diagnoses.age_at_diagnosis [days -> years], diagnoses
.ajcc_pathologic_tumor_stage) and cached as a sidecar TSV. The Cox protocol
is otherwise identical to tcga_multivariate_all.run() — constants and
parse_stage are imported from that module (no reimplementation):

- patient-level primary-tumour sample (matrix already one per patient;
  build_coad_gdc.py dedupes)
- OS from the build's data_clinical_patient_gdc.txt (GDC clinical API)
- signature (z) + AGE (z) + stage (ordinal 1-4, z); BH within 12 signatures

Output: data/tcga/COAD/COAD_cox_multivariate.json (+ sidecar clinical TSV)
"""
import importlib.util
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
TCGA = REPO / "data" / "tcga"
RAW = TCGA / "COAD" / "raw"
API = "https://api.gdc.cancer.gov"
SIDECAR = RAW / "data_clinical_patient_gdc_multiv.txt"

spec = importlib.util.spec_from_file_location(
    "tma", REPO / "scripts" / "tcga_multivariate_all.py")
tma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tma)  # module-level constants only; main() not run
SIGS = tma.SIGS
parse_stage = tma.parse_stage


def fetch_clinical():
    """Fetch per-case clinical from the GDC cases endpoint; cache sidecar.

    GDC Data Release 46 schema (verified 2026-09-27):
    - stage field is diagnoses.ajcc_pathologic_stage (the old
      ajcc_pathologic_tumor_stage name returns empty);
    - follow-up time lives in diagnoses.days_to_last_follow_up (the
      demographic-level field was retired), so OS for living patients uses
      the max diagnosis-level follow-up; dead patients use
      demographic.days_to_death;
    - age: diagnoses.age_at_diagnosis (days), fallback demographic
      .days_to_birth; both converted to years (/365.25).
    """
    if SIDECAR.exists():
        print("sidecar cached:", SIDECAR.name, flush=True)
        return pd.read_csv(SIDECAR, sep="\t", dtype=str)
    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({}))  # direct, no proxy
    rows, offset = {}, 0
    while True:
        q = {"filters": json.dumps(
                 {"op": "=", "content":
                  {"field": "project.project_id", "value": "TCGA-COAD"}}),
             "size": "100", "from": str(offset),
             "expand": "demographic,diagnoses",
             "fields": "submitter_id,vital_status",
             "format": "JSON"}
        req = urllib.request.Request(API + "/cases?" + urllib.parse.urlencode(q))
        with opener.open(req, timeout=120) as r:
            data = json.load(r)
        hits = data["data"]["hits"]
        if not hits:
            break
        for h in hits:
            dg = h.get("demographic") or {}
            dxs = h.get("diagnoses") or []
            stages = [d.get("ajcc_pathologic_stage")
                      for d in dxs if d.get("ajcc_pathologic_stage")]
            lfus = [d.get("days_to_last_follow_up")
                    for d in dxs if d.get("days_to_last_follow_up") is not None]
            aads = [d.get("age_at_diagnosis")
                    for d in dxs if d.get("age_at_diagnosis") is not None]
            vs = str(dg.get("vital_status") or "")
            dtd = dg.get("days_to_death")
            dtb = dg.get("days_to_birth")
            if vs.lower() == "dead" and dtd is not None:
                status, days = "1:DECEASED", float(dtd)
            elif vs.lower() == "alive":
                status = "0:LIVING"
                days = (float(max(lfus)) if lfus else None)
            else:
                status, days = None, None  # Not Reported
            age = None
            if aads:
                age = round(int(aads[0]) / 365.25, 2)
            elif dtb is not None:
                age = round(int(dtb) / -365.25, 2)
            rows[h["submitter_id"]] = {
                "PATIENT_ID": h["submitter_id"],
                "OS_STATUS": status or "",
                "OS_MONTHS": ("" if days is None
                              else round(days / 30.436875, 2)),
                "AGE": ("" if age is None else age),
                "AJCC_PATHOLOGIC_TUMOR_STAGE": stages[0] if stages else "",
            }
        offset += len(hits)
        if offset >= int(data["data"]["pagination"]["total"]):
            break
    df = pd.DataFrame(list(rows.values()))
    df.to_csv(SIDECAR, sep="\t", index=False)
    print(f"sidecar written: {len(df)} cases", flush=True)
    return df


def main():
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

    os_df = pd.read_csv(RAW / "data_clinical_patient_gdc.txt", sep="\t",
                        dtype=str, index_col=0)  # PATIENT_ID, OS_STATUS, OS_MONTHS
    side = fetch_clinical().set_index("PATIENT_ID")
    # DR46: the build-time demographic days_to_last_follow_up field was
    # retired, so the build's OS_MONTHS is empty for living patients; the
    # sidecar rebuilds OS from diagnosis-level follow-up. Cross-check events.
    cli = side[["OS_STATUS", "OS_MONTHS", "AGE",
                "AJCC_PATHOLOGIC_TUMOR_STAGE"]]
    ev_build = int(os_df["OS_STATUS"].astype(str).str.upper()
                   .str.contains("DECEASED").sum())
    ev_side = int(cli["OS_STATUS"].astype(str).str.upper()
                  .str.contains("DECEASED").sum())
    print(f"events cross-check: build clinical {ev_build} / sidecar {ev_side}",
          flush=True)
    if ev_build != ev_side:
        print("FATAL: event counts disagree — refusing to proceed", flush=True)
        return 1

    df = expr2.T.join(cli[["OS_STATUS", "OS_MONTHS", "AGE",
                           "AJCC_PATHOLOGIC_TUMOR_STAGE"]], how="inner")
    df["event"] = df["OS_STATUS"].astype(str).str.upper() \
        .str.contains("DECEASED").astype(int)
    df["dur"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
    df["age"] = pd.to_numeric(df["AGE"], errors="coerce")
    df["stage"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(parse_stage)
    df = df.dropna(subset=["dur", "age", "stage"])
    df = df[df["dur"] > 0]
    print(f"COAD: n={len(df)} events={int(df['event'].sum())}", flush=True)

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
    print(f"COAD multivariate: {n_sig}/{len(results)} BH-sig; HR {rng}",
          flush=True)
    out = {"cohort": "TCGA-COAD",
           "model": "multivariate Cox: signature (z) + AGE (z) + AJCC stage "
                    "(ordinal 1-4, z); BH within 12 signatures",
           "n": len(df), "events": int(df["event"].sum()),
           "source": ("GDC STAR matrix (build_coad_gdc.py) + GDC cases "
                      "endpoint clinical sidecar under Data Release 46 "
                      "schema (data_clinical_patient_gdc_multiv.txt; "
                      "ajcc_pathologic_stage / diagnosis-level "
                      "days_to_last_follow_up)"),
           "note": ("extends the SKCM multivariate sensitivity to this cancer "
                    "type; addresses the univariate-confounding attack on the "
                    "gradient"),
           "results": results}
    (TCGA / "COAD" / "COAD_cox_multivariate.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    print("written:", TCGA / "COAD" / "COAD_cox_multivariate.json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
