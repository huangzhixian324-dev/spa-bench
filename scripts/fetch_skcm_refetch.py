"""Refetch TCGA-SKCM for the sample-type sensitivity analysis — gene-limited
variant. The all-genes molecular-data request was rejected by the gateway
(502/503); the Cox pipeline only consumes the 12 immune signatures, so we
fetch exactly their genes (Entrez-resolved via POST /genes/fetch) for all
samples, plus patient-level OS_STATUS/OS_MONTHS.

Writes data/tcga/SKCM/refetch/{expression.tsv.gz, clinical.tsv}.

Faithfulness check: the all-tumour arm of the sensitivity pipeline must
reproduce the archived tcga_cox_all.json SKCM values (11/12 BH-significant;
HR 0.795-0.928) before the metastatic-only / primary-only arms are read.
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

BASE = "https://www.cbioportal.org/api"
H = {"Accept": "application/json", "Content-Type": "application/json"}
STUDY = "skcm_tcga_pan_can_atlas_2018"
OUT = Path(__file__).resolve().parents[1] / "data" / "tcga" / "SKCM" / "refetch"

SIG_GENES = sorted({g for s in [
    ['PDCD1','CD27','CTLA4','CD28','CD86','CD80','CD274','LAG3','HAVCR2','TIGIT','TNFRSF9','ICOS'],
    ['CD8A','CD8B'], ['CD19','CD79A','MS4A1'], ['NKG7','KLRD1','KLRF1'],
    ['IFNG','CXCL10','CXCL9','IDO1','STAT1'], ['GZMA','PRF1','GNLY'],
    ['PDCD1','CTLA4','LAG3','TIGIT','HAVCR2'], ['CXCL13','CCL19','CCL21'],
    ['FOXP3','IL2RA'], ['NOS2','IL12A','TNF'], ['CD14','CD33','ITGAM'],
    ['COL1A1','COL1A2','COL3A1','ACTA2','FAP']] for g in s})


def api_json(path, payload=None, timeout=900, retries=5):
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
            print(f"  attempt {attempt + 1}: {str(e)[:120]}", flush=True)
            time.sleep(15)
    raise RuntimeError(f"API failed: {path}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not (OUT / "expression.tsv.gz").exists():
        ent = {}
        for g in SIG_GENES:
            try:
                obj = api_json(f"/genes/{g}", timeout=60, retries=3)
                ent[g] = obj["entrezGeneId"]
            except Exception:
                print(f"  gene lookup failed: {g}", flush=True)
        missing = [g for g in SIG_GENES if g not in ent]
        print(f"entrez resolved: {len(ent)}/{len(SIG_GENES)}; "
              f"missing={missing}", flush=True)

        records = api_json(
            "/molecular-profiles/%s_rna_seq_v2_mrna/molecular-data/fetch"
            % STUDY,
            payload={"sampleListId": f"{STUDY}_all",
                     "entrezGeneIds": sorted(ent.values())})
        print(f"records = {len(records):,}", flush=True)

        ent2hugo = {v: k for k, v in ent.items()}
        mat = {}
        for rec in records:
            sym = ent2hugo.get(rec["entrezGeneId"])
            mat.setdefault(sym, {})[rec["sampleId"]] = rec.get("value")
        df = pd.DataFrame.from_dict(mat, orient="index")
        df.index.name = "Hugo_Symbol"
        df = df.reindex(sorted(df.columns), axis=1)
        df.to_csv(OUT / "expression.tsv.gz", sep="\t", compression="gzip")
        vals = df.values[~pd.isna(df.values)]
        print("expression: %s; value range %.3g..%.3g (median %.3g)" % (
            str(df.shape), vals.min(), vals.max(),
            float(pd.Series(vals).median())), flush=True)
    else:
        print("expression.tsv.gz already present; skipping fetch",
              flush=True)

    if not (OUT / "clinical.tsv").exists():
        rows = {}
        for dtype in ("SAMPLE", "PATIENT"):
            page, total = 0, 0
            while True:
                clin = api_json(
                    f"/studies/{STUDY}/clinical-data?clinicalDataType="
                    f"{dtype}&pageSize=20000&pageNumber={page}"
                    "&projection=DETAILED")
                if not clin:
                    break
                total += len(clin)
                for rec in clin:
                    attr = rec.get("clinicalAttributeId")
                    if attr in ("OS_STATUS", "OS_MONTHS"):
                        rows.setdefault(rec.get("patientId"), {})[attr] = \
                            rec.get("value")
                page += 1
            print(f"clinical[{dtype}] records: {total:,}", flush=True)
        cdf = pd.DataFrame.from_dict(rows, orient="index")
        cdf.index.name = "Patient Identifier"
        cdf.to_csv(OUT / "clinical.tsv", sep="\t")
        print("clinical:", cdf.shape, "OS_STATUS:",
              cdf["OS_STATUS"].value_counts().to_dict(), flush=True)
    else:
        print("clinical.tsv already present", flush=True)
    print("DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
