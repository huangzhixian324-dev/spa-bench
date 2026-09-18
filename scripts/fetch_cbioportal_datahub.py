"""Re-fetch the iAtlas expression matrices via the cBioPortal S3 DataHub
(the POST /api/studies/export endpoint used in 2026-09 was disabled, 405).

Downloads <study>.tar.gz from cbioportal-datahub.s3.amazonaws.com, extracts
the mRNA expression table, and writes the e1-compatible
expression.tsv.gz (index = Hugo_Symbol, columns = SAMPLE_ID) next to the
already-present clinical.tsv, matching the layout the e1 scripts expect.
"""
import io
import json
import sys
import tarfile
from pathlib import Path

import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "external"
STUDIES = ["mel_iatlas_liu_2019", "blca_iatlas_imvigor210_2017"]
BASE = "https://cbioportal-datahub.s3.amazonaws.com"

EXPR_CANDIDATES = [
    "data_mrna_seq_v2_rsem_Zscores.txt",
    "data_RNA_Seq_v2_mRNA_median_Zscores.txt",
    "data_mrna_seq_v2_rsem.txt",
    "data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt",
]


def fetch_one(study):
    url = f"{BASE}/{study}.tar.gz"
    print(f"[{study}] downloading {url}", flush=True)
    r = requests.get(url, timeout=900)
    r.raise_for_status()
    print(f"[{study}] {len(r.content) / 1e6:.1f} MB", flush=True)
    tf = tarfile.open(fileobj=io.BytesIO(r.content))
    names = tf.getnames()
    expr_name = None
    for cand in EXPR_CANDIDATES:
        hits = [n for n in names if n.endswith(cand)]
        if hits:
            expr_name = hits[0]
            break
    if expr_name is None:
        cand_expr = [n for n in names if "mrna" in n.lower()
                     or "rna_seq" in n.lower()]
        print(f"[{study}] no standard expression file; candidates: "
              f"{cand_expr}", flush=True)
        return None
    print(f"[{study}] extracting {expr_name}", flush=True)
    tf.extract(expr_name, OUT / study)
    src = OUT / study / expr_name
    df = pd.read_csv(src, sep="\t")
    # harmonize: Hugo_Symbol index, drop Entrez column, SAMPLE_ID columns
    if "Hugo_Symbol" in df.columns:
        df = df.drop(columns=[c for c in ("Entrez_Gene_Id",)
                              if c in df.columns])
        df = df.set_index("Hugo_Symbol")
    df = df[~df.index.duplicated(keep="first")]
    dst = OUT / study / "expression.tsv.gz"
    df.to_csv(dst, compression="gzip", sep="\t")
    print(f"[{study}] WROTE {dst}  shape={df.shape}  "
          f"first_cols={list(df.columns[:3])}", flush=True)
    return {"source": expr_name, "shape": list(df.shape)}


def main():
    summary = {}
    for study in STUDIES:
        try:
            summary[study] = fetch_one(study)
        except Exception as e:
            summary[study] = {"error": str(e)}
            print(f"[{study}] FAILED: {e}", flush=True)
    json.dump(summary, open(OUT / "expression_refetch_summary.json", "w"),
              indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
