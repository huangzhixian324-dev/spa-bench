"""Chunked re-fetch of the Liu 2019 / IMvigor210 expression matrices via the
cBioPortal molecular-data API. The single all-genes POST stalls server-side
(>1 h without response), so genes are fetched in chunks of 400 Entrez IDs
(~50 requests per study, seconds each) and assembled into the e1-compatible
expression.tsv.gz (index = Hugo_Symbol, columns = SAMPLE_ID).
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
CHUNK = 400


def api_json(path, payload=None, timeout=180, retries=5):
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
            time.sleep(10)
    raise RuntimeError(f"API failed: {path}")


def fetch_one(study):
    prof = f"{study}_rna_seq_mrna"
    # paginate the profile's gene list
    entrez, page = [], 0
    while True:
        genes = api_json(f"/molecular-profiles/{prof}/genes"
                         f"?pageSize=1000&pageNumber={page}")
        if not genes:
            break
        entrez += [g["entrezGeneId"] for g in genes]
        page += 1
    entrez = sorted(set(entrez))
    print(f"[{study}] profile={prof}, genes={len(entrez):,}", flush=True)

    rows, sample_order = {}, []
    for ci in range(0, len(entrez), CHUNK):
        chunk = entrez[ci:ci + CHUNK]
        recs = api_json(f"/molecular-profiles/{prof}/molecular-data/fetch",
                        payload={"sampleListId": f"{study}_all",
                                 "entrezGeneIds": chunk})
        for rec in recs:
            sym = rec.get("hugoGeneSymbol") or str(rec["entrezGeneId"])
            sid = rec["sampleId"]
            if sid not in sample_order:
                sample_order.append(sid)
            rows.setdefault(sym, {})[sid] = rec.get("value")
        print(f"    chunk {ci // CHUNK + 1}/{-(-len(entrez) // CHUNK)} "
              f"({len(recs):,} records)", flush=True)

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
    summary = {}
    for study in STUDIES:
        try:
            summary[study] = fetch_one(study)
        except Exception as e:
            summary[study] = {"error": str(e)}
            print(f"[{study}] FAILED: {e}", flush=True)
    json.dump(summary, open(OUT / "expression_chunkfetch_summary.json", "w"),
              indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
