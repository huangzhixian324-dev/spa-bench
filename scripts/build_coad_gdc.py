"""Rebuild TCGA-COAD for the Cox pipeline from the CURRENT GDC index.

The previously downloaded per-sample files are retired from the GDC index
(file UUIDs 404 / legacy 410), so this script re-downloads the current
TCGA-COAD primary-tumour STAR gene-expression quantifications (open
access) together with case/sample barcodes, and builds:

  data/tcga/COAD/raw/data_mrna_gdc_star_fpkm.txt   (Hugo_Symbol x barcode)
  data/tcga/COAD/raw/data_clinical_patient_gdc.txt (PATIENT_ID, OS_STATUS,
                                                    OS_MONTHS)
  data/tcga/COAD/raw/gdc_build_meta.json           (provenance)

Downloaded payloads are cached under data/tcga/TCGA-COAD_v2/ so a rerun
skips finished files.
"""
import concurrent.futures as cf
import json
import sys
from pathlib import Path

import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "data" / "tcga" / "TCGA-COAD_v2"
OUTD = REPO / "data" / "tcga" / "COAD" / "raw"
API = "https://api.gdc.cancer.gov"
DAYS_PER_MONTH = 30.436875
S = requests.Session()


def current_files():
    body = {
        "filters": {"op": "and", "content": [
            {"op": "=", "content":
             {"field": "cases.project.project_id", "value": "TCGA-COAD"}},
            {"op": "in", "content":
             {"field": "files.data_type",
              "value": ["Gene Expression Quantification"]}},
            {"op": "=", "content":
             {"field": "files.analysis.workflow_type",
              "value": "STAR - Counts"}},
            {"op": "in", "content":
             {"field": "cases.samples.sample_type",
              "value": ["Primary Tumor"]}},
        ]},
        "size": "1000",
        "expand": "cases.samples,cases.project",
        "fields": "file_id,file_name,file_size,cases.submitter_id,"
                  "cases.samples.submitter_id,cases.samples.sample_type",
        "format": "JSON",
    }
    r = S.post(f"{API}/files", data={k: (json.dumps(v) if isinstance(v, (dict, list)) else v)
                                     for k, v in body.items()}, timeout=120)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    print(f"current index: {len(hits)} primary-tumour STAR files")
    out = {}
    for h in hits:
        case = h["cases"][0]
        tumor = [s for s in case.get("samples", [])
                 if s["sample_type"] == "Primary Tumor"]
        if not tumor:
            continue
        bc = tumor[0]["submitter_id"]
        out[h["file_id"]] = {"barcode": bc,
                             "patient": case["submitter_id"],
                             "size": h.get("file_size", 0)}
    print(f"files with tumour barcode: {len(out)}")
    return out


def clinical_table():
    rows = {}
    offset = 0
    while True:
        r = S.get(f"{API}/cases", params={
            "filters": json.dumps({"op": "=", "content":
                                   {"field": "project.project_id",
                                    "value": "TCGA-COAD"}}),
            "size": "100", "from": str(offset), "expand": "demographic",
            "fields": "submitter_id,demographic.vital_status,"
                      "demographic.days_to_death,"
                      "demographic.days_to_last_follow_up",
            "format": "JSON"}, timeout=180)
        r.raise_for_status()
        hits = r.json()["data"]["hits"]
        if not hits:
            break
        for h in hits:
            pid = h["submitter_id"]
            dg = h.get("demographic") or {}
            vs = dg.get("vital_status")
            dtd = dg.get("days_to_death")
            dtlf = dg.get("days_to_last_follow_up")
            if vs is None:
                continue
            if str(vs).lower() == "dead":
                status, days = "1:DECEASED", dtd
            else:
                status, days = "0:LIVING", dtlf
            months = (None if days in (None, "", "NA")
                      else round(float(days) / DAYS_PER_MONTH, 2))
            rows[pid] = {"PATIENT_ID": pid, "OS_STATUS": status,
                         "OS_MONTHS": months}
        offset += len(hits)
        if offset >= int(r.json()["data"]["pagination"]["total"]):
            break
    print(f"clinical cases: {len(rows)}")
    return rows


def fetch(fid):
    dest = CACHE / f"{fid}.tsv"
    if dest.exists() and dest.stat().st_size > 10000:
        return fid, True
    for attempt in range(3):
        try:
            r = S.get(f"{API}/data/{fid}", timeout=300)
            r.raise_for_status()
            dest.write_bytes(r.content)
            # immediate structural validation: comment line + real header
            head = dest.read_text(errors="ignore")[:400]
            if "fpkm_unstranded" not in head:
                dest.unlink()
                raise ValueError("unexpected file structure")
            return fid, True
        except Exception as e:
            if attempt == 2:
                print(f"  [FAIL] {fid}: {e}", flush=True)
                return fid, False
    return fid, False


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    OUTD.mkdir(parents=True, exist_ok=True)
    files = current_files()
    # dedupe: one file per patient (smallest barcode first)
    per_patient = {}
    for fid, v in sorted(files.items(), key=lambda kv: kv[1]["barcode"]):
        per_patient.setdefault(v["patient"], (fid, v["barcode"]))
    print(f"patients to download: {len(per_patient)}")

    todo = [fid for fid, _ in per_patient.values()]
    S.headers.update({"Accept-Encoding": "gzip"})
    with cf.ThreadPoolExecutor(max_workers=12) as ex:
        for i, (fid, ok) in enumerate(ex.map(fetch, todo)):
            if (i + 1) % 50 == 0:
                print(f"  downloaded {i + 1}/{len(todo)}", flush=True)
    got = {fid for fid, ok in
           ((f, (CACHE / f"{f}.tsv").exists()) for f in todo) if ok}
    print(f"downloaded/cached: {len(got)}/{len(todo)}")

    # ---- expression matrix ----
    frames, used = [], {}
    for pid, (fid, bc) in sorted(per_patient.items()):
        if fid not in got:
            continue
        df = pd.read_csv(CACHE / f"{fid}.tsv", sep="\t", comment="#")
        df = df[df["gene_id"].astype(str).str.startswith("ENSG")]
        frames.append(pd.Series(df["fpkm_unstranded"].values,
                                index=df["gene_name"].astype(str).values,
                                name=bc))
        used[bc] = pid
    mat = pd.concat(frames, axis=1)
    mat = mat.groupby(level=0).mean()
    print(f"matrix: {mat.shape[0]} genes x {mat.shape[1]} samples")

    clin = clinical_table()
    keep = [bc for bc in mat.columns if used[bc] in clin]
    mat = mat[keep]
    clin_df = pd.DataFrame([clin[used[bc]] for bc in keep])
    print(f"with clinical: {mat.shape[1]} samples / "
          f"{int(clin_df['OS_STATUS'].str.contains('DECEASED').sum())} events")

    mat.to_csv(OUTD / "data_mrna_gdc_star_fpkm.txt", sep="\t",
               float_format="%.4f")
    clin_df.to_csv(OUTD / "data_clinical_patient_gdc.txt", sep="\t",
                   index=False)
    meta = {"n_files_current_index": len(files), "n_patients": len(used),
            "n_matrix_cols": int(mat.shape[1]),
            "events": int(clin_df["OS_STATUS"].str.contains("DECEASED")
                          .sum()),
            "source": "GDC current index: STAR - Counts, primary tumour, "
                      "fpkm_unstranded; clinical from GDC clinical API",
            "note": "COAD re-download replaces the withdrawn non-TCGA raw "
                    "files; SKCM/BRCA/LUAD unchanged (cBioPortal RSEM)"}
    json.dump(meta, open(OUTD / "gdc_build_meta.json", "w"), indent=1)
    print(json.dumps(meta, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
