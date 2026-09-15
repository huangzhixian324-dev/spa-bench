"""Gide cohort: splits generation (SHA-1 verified), h5ad export, and
fixed-scorer OOF AUROC fingerprints.

Fingerprints (benchmark_v33.json):
  GEP 0.830, TIDE 0.752, IMPRES 0.626, PD_L1 0.791,
  ElasticNet 0.630 (ML, seed/split dependent), ElasticNet_Var 0.605
  splits SHA-1[:16] = e0c8739e32e74797
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "cds_tool"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

API_DIR = REPO / "data/external/mel_iatlas_gide_2019_api"
OUT = REPO / "data/cohorts/Gide_2019_cBio/processed"
TARGET = "e0c8739e32e74797"


def bh_note():
    pass


def main():
    from sklearn.model_selection import StratifiedKFold, KFold

    clin = json.load(open(API_DIR / "clinical_data.json"))
    resp = {c["sampleId"]: c["value"] for c in clin
            if c["clinicalAttributeId"] == "RESPONSE"}
    compact = json.load(open(API_DIR / "expression_compact.json"))

    sym_of = {}
    df = pd.read_csv(REPO / "scripts/_hgnc_entrez.txt", sep="\t", dtype=str)
    for _, row in df.iterrows():
        sym_of[str(row["NCBI Gene ID(supplied by NCBI)"])] = \
            row["Approved symbol"]

    rows, symbols = [], []
    for entrez, rec in compact.items():
        sym = sym_of.get(entrez) or rec.get("hugo")
        if not sym:
            continue
        rows.append(rec["values"])
        symbols.append(sym)
    samples = sorted({s for r in rows for s in r})
    keep = [i for i, s in enumerate(samples)
            if s.endswith("_Pre") and s in resp and resp[s]]
    sample_ids = [samples[i] for i in keep]
    y = [1 if resp[s] in ("Complete Response", "Partial Response")
         else 0 for s in sample_ids]
    X_raw = np.array([[r.get(s, np.nan) for s in sample_ids] for r in rows]).T
    dfX = pd.DataFrame(X_raw, columns=symbols)
    dedup = dfX.T.groupby(level=0).max().T
    print(f"matrix: {dedup.shape}, responders {sum(y)}/{len(y)}")

    # ---- splits: try deterministic constructions against the registry hash
    hit = None
    y_arr = np.array(y)
    for skf_cls in (StratifiedKFold, KFold):
        for seed in (42, 0, 1, 123, 7, 2024, 2025, 2, 3, 4, 5, 10):
            try:
                sp = skf_cls(n_splits=5, shuffle=True, random_state=seed)
                folds = [{"fold": i, "train": t.tolist(), "test": e.tolist()}
                         for i, (t, e) in enumerate(
                             sp.split(np.zeros((len(y), 1)), y_arr))]
            except Exception:
                continue
            for label in ("sorted", "raw"):
                folds2 = [{"fold": f["fold"],
                           "train": sorted(f["train"]) if label == "sorted"
                           else f["train"],
                           "test": sorted(f["test"]) if label == "sorted"
                           else f["test"]} for f in folds]
                blob = json.dumps({"folds": folds2, "n": len(y), "k": 5})
                h = hashlib.sha1(blob.encode()).hexdigest()[:16]
                if h == TARGET:
                    hit = (folds2, f"{skf_cls.__name__} seed={seed} {label}")
                    print(f"*** SPLITS MATCH: {h} via {label} "
                          f"({skf_cls.__name__}, seed={seed})")
                    break
            if hit:
                break
        if hit:
            break
    if hit is None:
        print("no splits hash match across tried constructions; using "
              "StratifiedKFold(5, shuffle, seed=42) sorted as best guess")
        sp = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        folds = [{"fold": i, "train": sorted(t.tolist()),
                  "test": sorted(e.tolist())}
                 for i, (t, e) in enumerate(
                     sp.split(np.zeros((len(y), 1)), y_arr))]
    else:
        folds = hit[0]

    with open(OUT / "Gide_2019_cBio_splits.json", "w") as fh:
        fh.write(json.dumps({"folds": folds, "n": len(y), "k": 5}))

    # ---- h5ad export
    import scanpy as sc
    ad = sc.AnnData(X=dedup.values.astype(np.float64))
    ad.var_names = list(dedup.columns)
    ad.obs_names = [str(i) for i in range(len(y))]
    ad.obs["response"] = y
    ad.obs["sample_id"] = sample_ids
    ad.write_h5ad(OUT / "Gide_2019_cBio_processed.h5ad")
    print("wrote", OUT / "Gide_2019_cBio_processed.h5ad")

    # ---- fixed-scorer OOF fingerprints
    from base_method import get_method
    from sklearn.metrics import roc_auc_score
    bm = json.load(open(REPO / "results/benchmark/v33/"
                        "benchmark_v33.json"))["cohorts"]["Gide_2019_cBio"]
    results = {}
    for meth, key in [("IMPRES", "IMPRES"), ("GEP", "GEP"), ("TIDE", "TIDE"),
                      ("PD_L1_IHC", "PD_L1")]:
        try:
            m = get_method(meth)
            yt, yp = [], []
            for f in folds:
                preds, _ = m.fit_predict(ad, f["train"], f["test"])
                yt.extend(np.array(y)[f["test"]])
                yp.extend(preds)
            a = float(roc_auc_score(yt, yp))
            expect = bm[key]["auroc"]
            tag = ("*** MATCH" if abs(a - expect) < 0.001
                   else f"differs by {abs(a-expect):.4f}")
            results[meth] = {"rebuilt": round(a, 4),
                             "benchmark": round(expect, 4)}
            print(f"  {meth}: rebuilt {a:.4f} vs benchmark {expect:.4f}  {tag}")
        except Exception as e:
            print(f"  {meth}: FAILED {str(e)[:100]}")
    json.dump(results, open(OUT / "gide_fingerprint_check.json", "w"),
              indent=1)
    return 0


if __name__ == "__main__":
    main()
