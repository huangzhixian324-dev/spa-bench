"""Gene-chunked re-fetch of the Liu 2019 / IMvigor210 expression matrices.
The all-genes POST (any sample count) stalls server-side; gene-limited
requests work (validated on SKCM: 44 genes x 443 samples). Gene universe =
global /genes filtered to protein-coding (superset of any profile).
Chunks of 1,000 Entrez IDs; response records carry hugoGeneSymbol.
Writes the e1-compatible expression.tsv.gz per study.
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

BASE = "https://www.cbioportal.org/api"
H = {"Accept": "application/json", "Content-Type": "application/json"}
OUT = Path(__file__).resolve().parents[1] / "data" / "external"
STUDIES = ["mel_iatlas_liu_2019", "blca_iatlas_imvigor210_2017"]
CHUNK = 200


def api_json(path, payload=None, timeout=300, retries=5):
    for attempt in range(retries):
        try:
            if payload is None:
                req = urllib.request.Request(f"{BASE}{path}", headers=H)
            else:
                req = urllib.request.Request(
                    f"{BASE}{path}",
                    data=json.dumps(payload).encode(), headers=H)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"    attempt {attempt + 1}: {str(e)[:100]}", flush=True)
            time.sleep(8)
    raise RuntimeError(f"API failed: {path}")


def all_protein_coding():
    genes, page = [], 0
    while True:
        d = api_json(f"/genes?pageSize=10000&pageNumber={page}", timeout=120)
        if not d:
            break
        genes += [g["entrezGeneId"] for g in d
                  if g.get("type") == "protein-coding"]
        page += 1
    return sorted(set(genes))


def fetch_one(study, entrez):
    prof = f"{study}_rna_seq_mrna"
    rows = {}
    n_chunks = -(-len(entrez) // CHUNK)
    total_recs = 0
    for ci in range(0, len(entrez), CHUNK):
        chunk = entrez[ci:ci + CHUNK]
        recs = api_json(f"/molecular-profiles/{prof}/molecular-data/fetch",
                        payload={"sampleListId": f"{study}_all",
                                 "entrezGeneIds": chunk})
        total_recs += len(recs)
        for rec in recs:
            sym = rec.get("hugoGeneSymbol") or str(rec["entrezGeneId"])
            rows.setdefault(sym, {})[rec["sampleId"]] = rec.get("value")
        print(f"    gene-chunk {ci // CHUNK + 1}/{n_chunks} "
              f"({len(recs):,} records, cum {total_recs:,})", flush=True)

    df = pd.DataFrame.from_dict(rows, orient="index")
    df = df.reindex(columns=sorted(df.columns))
    df.index.name = "Hugo_Symbol"
    df = df[~df.index.isna()]
    df = df[~df.index.duplicated(keep="first")]
    dst = OUT / study / "expression.tsv.gz"
    df.to_csv(dst, sep="\t", compression="gzip")
    print(f"[{study}] WROTE {dst}  shape={df.shape}", flush=True)
    return {"profile": prof, "shape": [int(df.shape[0]), int(df.shape[1])]}


def main():
    entrez = all_protein_coding()
    print(f"protein-coding genes: {len(entrez):,}", flush=True)
    summary = {}
    for study in STUDIES:
        try:
            summary[study] = fetch_one(study, entrez)
        except Exception as e:
            summary[study] = {"error": str(e)}
            print(f"[{study}] FAILED: {e}", flush=True)
    json.dump(summary, open(OUT / "expression_chunkfetch_summary.json", "w"),
              indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
