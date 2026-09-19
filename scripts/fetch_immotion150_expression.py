"""Fetch the full IMmotion150 (iAtlas) TPM expression matrix via cBioPortal.

Profile: rcc_iatlas_immotion150_2018_rna_seq_mrna ("mRNA expression (TPM)")
Samples: all 263 (atezo arm 174 + sunitinib 89; arm split happens later)
Genes:   scripts/_hgnc_entrez.txt (43,687 Entrez-mapped HGNC symbols)

Batches of 1000 genes per POST /molecular-data/fetch; each batch is
parsed and written to raw/expr_batches/batch_NNNN.npz for resumability
(existing batches are skipped).  After all batches: merged into
raw/immotion150_tpm.npz {samples, genes, values}.

Usage: python scripts/fetch_immotion150_expression.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data/cohorts/IMmotion150/raw"
BATCH_DIR = RAW / "expr_batches"
PID = "rcc_iatlas_immotion150_2018_rna_seq_mrna"
SAMPLE_LIST = "rcc_iatlas_immotion150_2018_all"
BATCH = 1000

S = requests.Session()
S.headers["User-Agent"] = "spatbench-immotion150-fetch"


def main():
    BATCH_DIR.mkdir(parents=True, exist_ok=True)

    study = "rcc_iatlas_immotion150_2018"
    all_list = f"{study}_all"  # cBioPortal ALL-samples list convention (verified)
    sids = [x["sampleId"] for x in S.get(
        f"https://www.cbioportal.org/api/studies/{study}/samples"
        f"?projection=SUMMARY&page=0&pageSize=100000&direction=ASC",
        timeout=180).json()]
    sids = sorted(set(sids))
    print(f"samples: {len(sids)} (list {all_list})", flush=True)
    json.dump(sids, open(RAW / "samples.json", "w"))

    genes = []
    with open(REPO / "scripts/_hgnc_entrez.txt", encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2].isdigit():
                genes.append((parts[1], int(parts[2])))
    print(f"gene candidates: {len(genes)}", flush=True)

    n_batches = (len(genes) + BATCH - 1) // BATCH
    t0 = time.time()

    def do_batch(bi):
        chunk = genes[bi * BATCH:(bi + 1) * BATCH]
        out = BATCH_DIR / f"batch_{bi:04d}.npz"
        if out.exists():
            return bi, "cached"
        body = {"entrezGeneIds": [e for _, e in chunk],
                "sampleListId": SAMPLE_LIST}
        recs = None
        for attempt in range(5):
            try:
                r = S.post(
                    f"https://www.cbioportal.org/api/molecular-profiles/{PID}/molecular-data/fetch",
                    json=body, timeout=600)
                r.raise_for_status()
                recs = r.json()
                break
            except Exception as e:
                print(f"  batch {bi} attempt {attempt}: {str(e)[:70]}",
                      flush=True)
                time.sleep(3 + attempt * 5)
        if recs is None:
            return bi, "FAIL"
        by_gene = {}
        for x in recs:
            by_gene.setdefault(x["entrezGeneId"], {})[x["sampleId"]] = x["value"]
        e2i = {e: i for i, (_, e) in enumerate(chunk)}
        s2i = {s: i for i, s in enumerate(sids)}
        mat = np.full((len(chunk), len(sids)), np.nan, dtype=np.float32)
        for e, sv in by_gene.items():
            if e not in e2i:
                continue
            row = e2i[e]
            for s, v in sv.items():
                if s in s2i:
                    mat[row, s2i[s]] = v
        syms = [sym for sym, _ in chunk]
        eids = [e for _, e in chunk]
        np.savez_compressed(out, symbols=np.array(syms), entrez=np.array(eids),
                            values=mat)
        return bi, "ok"

    from concurrent.futures import ThreadPoolExecutor, as_completed
    done_ct = 0
    with ThreadPoolExecutor(4) as ex:
        futs = [ex.submit(do_batch, bi) for bi in range(n_batches)]
        for f in as_completed(futs):
            bi, status = f.result()
            done_ct += 1
            if status == "FAIL":
                print(f"  batch {bi} FAILED permanently", flush=True)
            if done_ct % 8 == 0 or done_ct == n_batches:
                el = time.time() - t0
                print(f"  progress {done_ct}/{n_batches} ({el/60:.1f} min)",
                      flush=True)

    # merge
    files = sorted(BATCH_DIR.glob("batch_*.npz"))
    print(f"merging {len(files)} batches ...", flush=True)
    syms, eids, mats = [], [], []
    for f in files:
        z = np.load(f, allow_pickle=True)
        syms.extend(z["symbols"].tolist())
        eids.extend(z["entrez"].tolist())
        mats.append(z["values"])
    M = np.concatenate(mats, axis=0)  # (G, S)
    np.savez_compressed(RAW / "immotion150_tpm.npz",
                        samples=np.array(sids), symbols=np.array(syms),
                        entrez=np.array(eids), values=M.T.astype(np.float32))
    n_present = int((~np.isnan(M)).any(axis=1).sum())
    print(f"MERGED: {M.shape[1]} samples x {M.shape[0]} gene-rows "
          f"({n_present} with data)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
