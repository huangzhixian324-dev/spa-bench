"""Fetch and persist the full Gide 2019 expression matrix (all genes x all
samples) from the cBioPortal API in one call, streaming to disk.
"""
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

BASE = "https://www.cbioportal.org/api"
H = {"Accept": "application/json", "Content-Type": "application/json"}
STUDY = "mel_iatlas_gide_2019"
PROFILE = "mel_iatlas_gide_2019_rna_seq_mrna"
OUT = Path(__file__).resolve().parents[1] / "data/external/mel_iatlas_gide_2019_api"
OUT.mkdir(parents=True, exist_ok=True)
RAW_JSON = OUT / "molecular_data_all.json"


def main():
    for attempt in range(5):
        try:
            req = urllib.request.Request(
                f"{BASE}/molecular-profiles/{PROFILE}/molecular-data/fetch",
                data=json.dumps({"sampleListId": f"{STUDY}_all"}).encode(),
                headers=H)
            with urllib.request.urlopen(req, timeout=1800) as r:
                blob = r.read()
            print(f"fetched {len(blob)/1e6:.1f} MB json", flush=True)
            break
        except Exception as e:
            print(f"attempt {attempt}: {str(e)[:80]}", flush=True)
            time.sleep(5)
    else:
        return 1
    tmp = str(RAW_JSON) + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(blob)
    os.replace(tmp, RAW_JSON)
    print("saved", RAW_JSON, f"{RAW_JSON.stat().st_size/1e6:.1f} MB", flush=True)

    # parse into a compact matrix file (gene -> {sample: value})
    records = json.load(open(RAW_JSON, encoding="utf-8"))
    ent2hugo, ent2sym = {}, {}
    mat = {}
    for rec in records:
        e = rec["entrezGeneId"]
        h = ent2hugo.setdefault(e, rec.get("hugoGeneSymbol"))
        mat.setdefault(e, {})[rec["sampleId"]] = rec.get("value")
    print(f"genes: {len(mat)}", flush=True)
    compact = {str(e): {"hugo": h, "values": v}
               for e, (h, v) in ((e, (ent2hugo[e], mat[e])) for e in mat)}
    json.dump(compact, open(OUT / "expression_compact.json", "w"))
    print("saved expression_compact.json", flush=True)
    os.remove(RAW_JSON)
    print("DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
