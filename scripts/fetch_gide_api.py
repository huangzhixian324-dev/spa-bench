"""Re-fetch the Gide 2019 (mel_iatlas_gide_2019) expression + clinical data
from the cBioPortal API (the /studies/export endpoint is disabled on the
public instance, so we page through molecular-data instead).

Outputs (data/external/mel_iatlas_gide_2019_api/):
  samples.json        - all samples with patientId / sampleId
  clinical_data.json  - per-sample clinical attributes (RECIST etc.)
  expression.jsonl    - one line per gene: {entrez, hugo, {sample: value}}
Checkpointed per batch so an interrupted run can resume.
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

BASE = "https://www.cbioportal.org/api"
HEADERS = {"Accept": "application/json", "Content-Type": "application/json"}
STUDY = "mel_iatlas_gide_2019"
PROFILE = "mel_iatlas_gide_2019_rna_seq_mrna"
OUT = Path(__file__).resolve().parents[1] / "data/external/mel_iatlas_gide_2019_api"
OUT.mkdir(parents=True, exist_ok=True)


def api(path, params=None, body=None, retries=5):
    url = f"{BASE}{path}"
    data = json.dumps(body).encode() if body is not None else None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, data=data,
                                         headers=HEADERS)
            sep = "&" if params and "?" in url else "?"
            if params:
                from urllib.parse import urlencode
                url_full = url + sep + urlencode(params)
            else:
                url_full = url
            req = urllib.request.Request(url_full, data=data, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            print(f"  retry {i}: {str(e)[:70]}", flush=True)
            time.sleep(2 + i * 3)
    raise RuntimeError(f"api failed: {path}")


def main():
    samples = api(f"/studies/{STUDY}/samples",
                  {"pageSize": 100000, "pageNumber": 0})
    sample_ids = [s["sampleId"] for s in samples]
    print(f"samples: {len(sample_ids)}", flush=True)
    json.dump(samples, open(OUT / "samples.json", "w"), indent=0)

    clin = api(f"/studies/{STUDY}/clinical-data",
               {"clinicalDataType": "SAMPLE", "projection": "DETAILED",
                "pageSize": 100000, "pageNumber": 0})
    json.dump(clin, open(OUT / "clinical_data.json", "w"), indent=0)
    recist = [c for c in clin
              if c.get("clinicalAttributeId", "").upper().find("RECIST") >= 0
              or c.get("clinicalAttributeId", "").upper() == "RESPONSE"]
    print(f"clinical attrs: {len(clin)}; RECIST-like: {len(recist)}", flush=True)

    # gene list of the profile
    genes = []
    page = 0
    while True:
        batch = api(f"/molecular-profiles/{PROFILE}/genes",
                    {"pageSize": 1000, "pageNumber": page,
                     "projection": "SUMMARY"})
        genes.extend(batch)
        if len(batch) < 1000:
            break
        page += 1
    entrez = [g["entrezGeneId"] for g in genes if g.get("entrezGeneId")]
    print(f"profile genes: {len(genes)} (entrez: {len(entrez)})", flush=True)
    json.dump(genes, open(OUT / "profile_genes.json", "w"), indent=0)

    # expression in batches of 400 genes x all samples
    done = set()
    jsonl = OUT / "expression.jsonl"
    if jsonl.exists():
        for line in jsonl.read_text(encoding="utf-8").splitlines():
            try:
                done.add(json.loads(line)["entrezGeneId"])
            except Exception:
                pass
    with open(jsonl, "a", encoding="utf-8") as fh:
        for i in range(0, len(entrez), 400):
            batch = entrez[i:i + 400]
            if any(e in done for e in batch):
                continue
            res = api(f"/molecular-profiles/{PROFILE}/molecular-data/fetch",
                      body={"sampleListId": f"{STUDY}_all",
                            "entrezGeneIds": batch})
            fh.write(json.dumps(res) + "\n")
            fh.flush()
            print(f"  expression batch {i // 400 + 1}/"
                  f"{(len(entrez) + 399) // 400} done", flush=True)
    print("DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
