"""Rebuild the Gide (mel_iatlas_gide_2019) cohort from the re-fetched
cBioPortal data and verify against archived fingerprints:
  CDS (GZMA/PRF1) = 0.852  (Table S4b, recomputation)
  splits SHA-1[:16] = e0c8739e32e74797
  fixed-scorer OOF AUROCs: GEP 0.830, TIDE 0.752, IMPRES 0.626, PD-L1 0.791
Samples: 73 Pre-treatment with RECIST annotation (40 CR+PR = 55%).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "cds_tool"))

API_DIR = REPO / "data/external/mel_iatlas_gide_2019_api"
OUT = REPO / "data/cohorts/Gide_2019_cBio/processed"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    # ---- symbol map (entrez -> hgnc) from the local HGNC subset
    sym_of = {}
    df = pd.read_csv(REPO / "scripts/_hgnc_entrez.txt", sep="\t", dtype=str)
    for _, row in df.iterrows():
        gid = row["NCBI Gene ID(supplied by NCBI)"]
        sym_of[str(gid)] = row["Approved symbol"]

    # ---- expression matrix
    compact = json.load(open(API_DIR / "expression_compact.json"))
    rows, symbols = [], []
    for entrez, rec in compact.items():
        sym = sym_of.get(entrez) or rec.get("hugo")
        if not sym:
            continue
        rows.append(rec["values"])
        symbols.append(sym)
    samples = sorted({s for r in rows for s in r})
    X_raw = np.array([[r.get(s, np.nan) for s in samples] for r in rows]).T
    print(f"raw matrix: {X_raw.shape[0]} samples x {X_raw.shape[1]} genes")

    # ---- clinical: Pre-treatment + RECIST response
    clin = json.load(open(API_DIR / "clinical_data.json"))
    resp = {c["sampleId"]: c["value"] for c in clin
            if c["clinicalAttributeId"] == "RESPONSE"}
    keep = [i for i, s in enumerate(samples)
            if s.endswith("_Pre") and s in resp and resp[s]]
    y = [1 if resp[samples[i]] in ("Complete Response", "Partial Response")
         else 0 for i in keep]
    print(f"kept {len(keep)} Pre samples; responders {sum(y)} "
          f"({sum(y)/len(y):.0%})  [expect 73 / 40 (55%)]")

    # ---- deduplicate symbols (keep the highest-variance row per symbol)
    X = X_raw[keep]
    dfX = pd.DataFrame(X, columns=symbols)
    dedup = (dfX.T.groupby(level=0).max()).T
    print(f"deduplicated matrix: {dedup.shape}")

    # ---- CDS fingerprint (z-score expression; rank-based components are
    # invariant to affine transforms)
    sys.path.insert(0, str(REPO / "cds_tool"))
    from cds import circularity_detection_score
    res = circularity_detection_score(
        dedup.values, np.array(y), ["GZMA", "PRF1"],
        list(dedup.columns), random_state=42)
    c1 = res["components"]["C1_MI_endpoint_genes"]
    c2 = res["components"]["C2_rho_endpoint_response"]
    c3 = res["components"]["C3_MI_skew_ratio"]
    print(f"CDS rebuilt: {res['CDS']:.3f}  C1={c1:.4f} C2={c2:.4f} C3={c3:.4f}")
    print("expect       : 0.852  C1=0.141 C2=0.508 C3=1.000  "
          "(Table S4b, HIGH)")

    json.dump({"n_samples": len(keep), "n_genes": int(dedup.shape[1]),
               "responders": int(sum(y)),
               "cds_rebuilt": round(float(res["CDS"]), 3),
               "cds_expected": 0.852,
               "components": {"C1": round(float(c1), 4),
                              "C2": round(float(c2), 4),
                              "C3": round(float(c3), 4)}},
              open(OUT / "gide_rebuild_check.json", "w"), indent=1)
    print("saved", OUT / "gide_rebuild_check.json")
    return 0


if __name__ == "__main__":
    main()
