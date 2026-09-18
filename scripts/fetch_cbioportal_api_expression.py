"""Re-fetch the Liu 2019 / IMvigor210 expression matrices via the cBioPortal
molecular-data API (the study-export endpoint used previously is disabled,
405; the S3 datahub excludes iAtlas studies). Mirrors scripts/
fetch_gide_full.py, which produced the archived Gide matrix the same way.

For each study: probe the molecular profile id, then POST
/api/molecular-profiles/{profile}/molecular-data/fetch with the all-sample
list, and write the e1-compatible expression.tsv.gz
(index = Hugo_Symbol, columns = SAMPLE_ID) under data/external/<study>/.
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


def api_json(path, payload=None, timeout=1800, retries=4):
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
            print(f"  attempt {attempt + 1}: {str(e)[:100]}", flush=True)
            time.sleep(8)
    raise RuntimeError(f"API failed: {path}")


def fetch_one(study):
    profiles = api_json(f"/studies/{study}/molecular-profiles")
    mrna = [p for p in profiles
            if "mrna" in p["molecularProfileId"].lower()
            or "rna_seq" in p["molecularProfileId"].lower()]
    if not mrna:
        print(f"[{study}] no mRNA profile: "
              f"{[p['molecularProfileId'] for p in profiles]}", flush=True)
        return {"error": "no mRNA profile"}
    prof = mrna[0]["molecularProfileId"]
    print(f"[{study}] profile = {prof}", flush=True)

    records = api_json(f"/molecular-profiles/{prof}/molecular-data/fetch",
                       payload={"sampleListId": f"{study}_all"})
    print(f"[{study}] records = {len(records):,}", flush=True)
    mat, hugo = {}, {}
    for rec in records:
        e = rec["entrezGeneId"]
        hugo.setdefault(e, rec.get("hugoGeneSymbol"))
        mat.setdefault(e, {})[rec["sampleId"]] = rec.get("value")
    df = pd.DataFrame.from_dict(mat, orient="index")
    df.index = [hugo[e] for e in df.index]
    df = df[~df.index.isna()]
    df.index = df.index.astype(str)
    df = df[~df.index.duplicated(keep="first")]
    dst = OUT / study / "expression.tsv.gz"
    df.to_csv(dst, compression="gzip", sep="\t")
    print(f"[{study}] WROTE {dst}  shape={df.shape}  "
          f"first_cols={list(df.columns[:3])}", flush=True)
    return {"profile": prof, "shape": [int(df.shape[0]), int(df.shape[1])]}


def main():
    summary = {}
    for study in STUDIES:
        try:
            summary[study] = fetch_one(study)
        except Exception as e:
            summary[study] = {"error": str(e)}
            print(f"[{study}] FAILED: {e}", flush=True)
    p = OUT / "expression_apifetch_summary.json"
    json.dump(summary, open(p, "w"), indent=1)
    print("summary ->", p, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
