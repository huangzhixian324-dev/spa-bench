"""Build the IMmotion150 (iAtlas) ccRCC cohort for the executed E1.

Data sources:
  - Expression: cBioPortal API molecular-data/fetch, profile
    rcc_iatlas_immotion150_2018_rna_seq_mrna ("mRNA expression (TPM)"),
    263 samples x ~40k Entrez genes, fetched in 1000-gene batches into
    raw/expr_batches/*.npz and merged to raw/immotion150_tpm.npz by
    scripts/fetch_immotion150_expression.py.
  - Clinical: same iAtlas study, per-patient /clinical-data endpoint
    (263/263 fetched; snapshot in data/external/_immotion150_clinical_probe.json).
    Fields: RESPONDER (RECIST CR/PR = TRUE), PFS_STATUS, PFS_MONTHS,
    ICI_RX (Atezolizumab 174 / None 89).

Cohort definition:
  - Atezolizumab arm only (ICI_RX = 'Atezolizumab'), n = 174;
    9 patients lack RESPONDER and are dropped from the RECIST endpoint
    (retained for PFS endpoints).
  - Note: iAtlas does not separate atezo monotherapy vs atezo+bevacizumab
    randomised arms (SAMPLE_TREATMENT empty) -- the 174-patient ICI arm
    pools both; disclosed as a caveat.

Endpoints:
  - RECIST: RESPONDER TRUE=1 / FALSE=0 (CR/PR vs SD/PD, the trial's own
    RECIST coding; identical to the SPATBench convention).
  - DCB6mo (primary PFS): PFS >= 6 months AND PFS_STATUS = not progressed
    = 1 (durable benefit; the manuscript's Jung-2019 convention -- this
    cohort HAS an event indicator, unlike GSE274975).
  - PFS6mo (sensitivity, GSE274975-comparable): PFS < 6 months = early
    progressor, time-only.

Processing:
  - TPM -> log2(TPM+1); duplicate HGNC symbols aggregated by MEAN
    (TPM is not additive, unlike raw counts).
  - Splits: stratified k = min(5, n_pos, n_neg), seed 42, SHA-1 recorded.

Output:
  data/cohorts/IMmotion150/processed/IMmotion150_HGNC.h5ad
  data/cohorts/IMmotion150/processed/IMmotion150_splits.json
  data/cohorts/IMmotion150/processed/build_qc.json
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data/cohorts/IMmotion150/raw"
PROC = REPO / "data/cohorts/IMmotion150/processed"
CLIN = REPO / "data/external/_immotion150_clinical_probe.json"
SEED = 42


def stratified_folds(y, k, seed):
    rng = np.random.RandomState(seed)
    fold_of = np.full(len(y), -1)
    for cls in np.unique(y):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for pos, i in enumerate(idx):
            fold_of[i] = pos % k
    return [{"train": [int(i) for i in np.where(fold_of != f)[0]],
             "test": [int(i) for i in np.where(fold_of == f)[0]]}
            for f in range(k)]


def main():
    PROC.mkdir(parents=True, exist_ok=True)

    z = np.load(RAW / "immotion150_tpm.npz", allow_pickle=True)
    samples = z["samples"].tolist()
    symbols = z["symbols"].tolist()
    entrez = z["entrez"].tolist()
    values = z["values"]  # (S, G) TPM, NaN = missing
    print(f"raw matrix: {values.shape[0]} samples x {values.shape[1]} genes")

    clin = json.load(open(CLIN))

    # ---- atezo arm ----
    keep, resp, pfs_m, pfs_e = [], [], [], []
    for s in samples:
        c = clin.get(s, {})
        if c.get("ICI_RX") != "Atezolizumab":
            continue
        keep.append(s)
        r = c.get("RESPONDER")
        resp.append(1 if r == "TRUE" else (0 if r == "FALSE" else -1))
        pm = c.get("PFS_MONTHS")
        pfs_m.append(float(pm) if pm not in (None, "", "NA") else np.nan)
        ps = c.get("PFS_STATUS", "")
        pfs_e.append(1 if str(ps).startswith("1") else (0 if str(ps).startswith("0") else -1))
    keep = np.array(keep)
    resp = np.array(resp, dtype=int)
    pfs_m = np.array(pfs_m, dtype=float)
    pfs_e = np.array(pfs_e, dtype=int)
    print(f"atezo arm: {len(keep)} patients; RESPONDER "
          f"TRUE={int((resp==1).sum())} FALSE={int((resp==0).sum())} "
          f"missing={int((resp<0).sum())}")

    X = values[[samples.index(s) for s in keep], :]
    n_nan_samples = int(np.isnan(X).all(axis=1).sum())
    print(f"expression subset: {X.shape}; all-NaN samples: {n_nan_samples}")

    # ---- gene aggregation: duplicate symbols -> MEAN (TPM not additive) ----
    ok = ~np.isnan(X).all(axis=0)          # genes with any data
    syms = np.array(symbols)[ok]
    Xg = X[:, ok]
    uniq, inv = np.unique(syms, return_inverse=True)
    G = len(uniq)
    agg = np.zeros((Xg.shape[0], G), dtype=np.float64)
    cnt = np.zeros((Xg.shape[0], G), dtype=np.float32)
    for j in range(Xg.shape[1]):
        g = inv[j]
        col = Xg[:, j]
        m = ~np.isnan(col)
        agg[m, g] += col[m]
        cnt[m, g] += 1
    with np.errstate(invalid="ignore"):
        agg = agg / np.maximum(cnt, 1)
    print(f"unique symbols: {G} (dup rows {Xg.shape[1] - G} aggregated by mean)")

    # ---- log2(TPM+1) ----
    agg = np.log2(agg + 1.0)
    Xp = agg.astype(np.float32)

    # ---- endpoints ----
    recist_ok = resp >= 0
    response = resp.copy()
    response[~recist_ok] = 0  # placeholder; RECIST analyses use mask
    n_res = int((resp == 1).sum())
    n_non = int((resp == 0).sum())
    dcb6 = ((pfs_m >= 6) & (pfs_e == 0)).astype(int)
    pfs6 = (pfs_m < 6).astype(int)
    print(f"RECIST: responders {n_res}/{n_res + n_non}")
    print(f"DCB6mo (PFS>=6 & no event): {int(dcb6.sum())}/{len(dcb6)}")
    print(f"PFS<6mo (time-only): {int(pfs6.sum())}/{len(pfs6)}")

    # ---- AnnData ----
    import anndata as ad
    obs = ad.AnnData(X=Xp).obs
    obs = obs.set_index(np.array([f"{s}" for s in keep], dtype=object))
    df_obs = {
        "response": response,
        "recist_available": recist_ok.astype(int),
        "responder_true": (resp == 1).astype(int),
        "responder_false": (resp == 0).astype(int),
        "pfs_months": pfs_m,
        "pfs_event": pfs_e,
        "dcb6": dcb6,
        "pfs6_timeonly": pfs6,
        "ici_rx": ["Atezolizumab"] * len(keep),
        "cancer_type": ["ccRCC (IMmotion150)"] * len(keep),
        "patient_id": keep.astype(object),
    }
    import pandas as pd
    obs_df = pd.DataFrame(df_obs, index=[str(s) for s in keep])
    adata = ad.AnnData(X=Xp, obs=obs_df)
    adata.var_names = [str(u) for u in uniq]
    adata.uns["provenance"] = (
        "iAtlas harmonized IMmotion150 via cBioPortal API "
        "(rcc_iatlas_immotion150_2018_rna_seq_mrna, TPM; clinical per-patient "
        "endpoint). atezo arm pooled (mono+bev arms not separable in iAtlas). "
        "log2(TPM+1); duplicate symbols mean-aggregated")

    # ---- splits (RECIST labels) ----
    y = resp[recist_ok].astype(int)
    idx_ok = np.where(recist_ok)[0]
    k = min(5, int((y == 1).sum()), int((y == 0).sum()))
    folds = stratified_folds(y, k, SEED)
    # remap fold indices to the full-cohort rows
    folds_full = [{"train": [int(idx_ok[i]) for i in f["train"]],
                   "test": [int(idx_ok[i]) for i in f["test"]]} for f in folds]
    payload = {"cohort": "IMmotion150", "k": k, "seed": SEED,
               "endpoint": f"RECIST (responders {n_res}/{n_res + n_non})",
               "folds": folds_full,
               "recist_available_rows": [int(i) for i in idx_ok]}
    payload["splits_sha1_16"] = hashlib.sha1(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]

    adata.write_h5ad(PROC / "IMmotion150_HGNC.h5ad")
    json.dump(payload, open(PROC / "IMmotion150_splits.json", "w"), indent=1)
    qc = {"n_patients_atezo": int(len(keep)),
          "n_recist_ok": int(recist_ok.sum()), "n_responders": n_res,
          "n_dcb6": int(dcb6.sum()), "n_pfs6_timeonly": int(pfs6.sum()),
          "n_unique_symbols": int(G), "splits_sha1_16": payload["splits_sha1_16"]}
    json.dump(qc, open(PROC / "build_qc.json", "w"), indent=1)
    print("\nQC:", json.dumps(qc, indent=1))
    print(f"\nWROTE {PROC / 'IMmotion150_HGNC.h5ad'}")
    print(f"WROTE {PROC / 'IMmotion150_splits.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
