"""Rebuild the Jung (GSE135222) cohort: Ensembl -> HGNC via Ensembl REST,
DCB labels (PFS >= 180 d AND no progression event = 6/27), splits with
SHA-1 verification against run_registry, then fixed-scorer OOF AUROC
fingerprints against benchmark_v33.json.

Expected fingerprints (benchmark_v33.json / manuscript Table 6):
  IMPRES 0.583, GEP 0.786, TIDE 0.278, PD_L1 0.659
  splits SHA-1[:16] = 9cd158fbd9e77626 (HGNC)
"""
import gzip
import hashlib
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
os.chdir(REPO)
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

RAW = REPO / "data/cohorts/Jung_2019/raw"
TARGET_SHAS = {"9cd158fbd9e77626", "6821301b46ccfb99"}
HEADERS = {"Content-Type": "application/json", "Accept": "application/json",
           "User-Agent": "spatbench-rebuild@example.org"}


def ensembl_map(ids):
    out = {}
    for i in range(0, len(ids), 900):
        batch = ids[i:i + 900]
        body = json.dumps({"ids": batch}).encode()
        for attempt in range(5):
            try:
                req = urllib.request.Request(
                    "https://rest.ensembl.org/lookup/id", data=body,
                    headers=HEADERS)
                with urllib.request.urlopen(req, timeout=60) as r:
                    res = json.load(r)
                for e, info in res.items():
                    if info and info.get("display_name"):
                        out[e] = info["display_name"]
                break
            except Exception as e:
                print(f"  batch {i} retry {attempt}: {str(e)[:60]}")
                time.sleep(2 + attempt * 3)
        print(f"  mapped {len(out)}/{len(ids)}")
    return out


def load_raw():
    genes, samples, mat = [], [], []
    with gzip.open(RAW / "GSE135222_exp.tsv.gz", "rt") as f:
        header = f.readline().rstrip().split("\t")[1:]
        samples = [s.strip() for s in header]
        for line in f:
            parts = line.rstrip().split("\t")
            genes.append(parts[0])
            mat.append([float(x) for x in parts[1:]])
    titles, flags, times = [], [], []
    with gzip.open(RAW / "GSE135222_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"').replace(" ", "")
                          for t in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and \
                    "progression-free survival (pfs):" in line:
                flags = [t.strip('"').split(": ")[1]
                         for t in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and \
                    '\t"pfs.time' in line:
                times = [t.strip('"').split(": ")[1]
                         for t in line.rstrip().split("\t")[1:]]
    title2idx = {t: i for i, t in enumerate(titles)}
    # order expression columns to match series-matrix title order
    col_order = [samples.index(t) for t in titles]
    y = [1 if (int(times[title2idx[t]]) >= 180 and
               flags[title2idx[t]] == "0") else 0 for t in titles]
    return genes, mat, col_order, titles, y


def main():
    genes, mat, col_order, titles, y = load_raw()
    y_arr = np.array(y)
    print(f"n={len(titles)} DCB={int(y_arr.sum())} (expect 6)")
    assert y_arr.sum() == 6, "DCB label mismatch"

    mat_arr = np.array(mat)[:, col_order].T  # (27, n_genes) TPM
    print("log2(TPM+1)...")

    cache = REPO / "scripts/_jung_ensembl_map.json"
    if cache.exists():
        emap = json.load(open(cache))
    else:
        emap = {}
    if len(emap) < len(set(g.split(".")[0] for g in genes)):
        print(f"Ensembl REST mapping for {len(genes)} ids ...", flush=True)
        clean = [g.split(".")[0] for g in genes]
        uniq = sorted(set(clean) - set(emap))
        for i in range(0, len(uniq), 900):
            batch = uniq[i:i + 900]
            body = json.dumps({"ids": batch}).encode()
            for attempt in range(5):
                try:
                    req = urllib.request.Request(
                        "https://rest.ensembl.org/lookup/id", data=body,
                        headers=HEADERS)
                    with urllib.request.urlopen(req, timeout=120) as r:
                        res = json.load(r)
                    for e, info in res.items():
                        if info and info.get("display_name"):
                            emap[e] = info["display_name"]
                    json.dump(emap, open(cache, "w"), indent=0)
                    print(f"  batch {i // 900 + 1}: total mapped "
                          f"{len(emap)}/{len(clean)}", flush=True)
                    break
                except Exception as e:
                    print(f"  batch {i // 900 + 1} retry {attempt}: "
                          f"{str(e)[:60]}", flush=True)
                    time.sleep(2 + attempt * 3)
    print(f"mapped {len(emap)}/{len(genes)} ensembl ids to symbols", flush=True)

    clean = [g.split(".")[0] for g in genes]
    keep = [i for i, e in enumerate(clean) if e in emap]
    symbols = [emap[clean[i]] for i in keep]
    order = np.argsort(symbols)
    symbols_sorted = [symbols[i] for i in order]
    X = np.log2(mat_arr[:, keep][:, order] + 1.0)
    print(f"matrix: {X.shape[0]} x {X.shape[1]} HGNC genes")

    from sklearn.model_selection import StratifiedKFold, KFold
    hit = None
    for seed, skf_cls in itertools_product():
        try:
            splitter = skf_cls(n_splits=5, shuffle=True, random_state=seed)
            folds = [{"fold": i, "train": t.tolist(), "test": e.tolist()}
                     for i, (t, e) in enumerate(
                         splitter.split(np.zeros((len(y), 1)), y_arr))]
            blob = json.dumps({"folds": folds, "n": len(y), "k": 5})
            h = hashlib.sha1(blob.encode()).hexdigest()[:16]
        except Exception:
            continue
        print(f"  {skf_cls.__name__} seed={seed}: {h}")
        if h in TARGET_SHAS:
            hit = (seed, blob, folds)
            print(f"  *** MATCH {h} (seed={seed})")
            break
    if hit is None:
        print("no split hash match; still writing matrix + best-guess splits")
        seed, blob, folds = 42, None, None
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        folds = [{"fold": i, "train": t.tolist(), "test": e.tolist()}
                 for i, (t, e) in enumerate(
                     skf.split(np.zeros((len(y), 1)), y_arr))]

    out = REPO / "data/cohorts/Jung_2019/processed"
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "Jung_2019_HGNC_X.npy", X)
    json.dump({"symbols": symbols_sorted, "obs": titles, "response": y},
              open(out / "Jung_2019_HGNC_meta.json", "w"), indent=0)
    with open(out / "Jung_2019_HGNC_splits.json", "w") as fh:
        fh.write(json.dumps({"folds": folds, "n": len(y), "k": 5}))
    print("wrote matrix/meta/splits under", out)

    # fixed-scorer fingerprints
    try:
        from base_method import get_method  # noqa
        import scanpy as sc
        ad = sc.AnnData(X=X)
        ad.var_names = symbols_sorted
        ad.obs["response"] = y
        bm = json.load(open(REPO / "results/benchmark/v33/"
                            "benchmark_v33.json"))["cohorts"]["Jung_2019_DCB"]
        for meth, expect in [("IMPRES", bm["IMPRES"]["auroc"]),
                             ("GEP", bm["GEP"]["auroc"]),
                             ("PD_L1_IHC", bm["PD_L1"]["auroc"])]:
            m = get_method(meth)
            yt, yp = [], []
            for f in folds:
                preds, _ = m.fit_predict(ad, f["train"], f["test"])
                yt.extend(np.array(y)[f["test"]])
                yp.extend(preds)
            from sklearn.metrics import roc_auc_score
            a = roc_auc_score(yt, yp)
            print(f"  {meth}: rebuilt {a:.4f} vs benchmark {expect:.4f} "
                  f"{'*** MATCH' if abs(a - expect) < 0.001 else 'differs'}")
    except Exception as e:
        print("fingerprint run failed:", str(e)[:120])
    return 0


def itertools_product():
    from sklearn.model_selection import StratifiedKFold, KFold
    for seed in [42, 0, 1, 123, 7, 2024]:
        yield seed, StratifiedKFold
    for seed in [42, 0, 1, 123]:
        yield seed, KFold


if __name__ == "__main__":
    raise SystemExit(main())
