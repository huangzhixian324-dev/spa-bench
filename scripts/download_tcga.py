#!/usr/bin/env python3
"""Download TCGA expression + survival data for SPA-Bench."""
import os, json, time, requests, urllib3, gzip, pandas as pd
from pathlib import Path
urllib3.disable_warnings()
os.environ['http_proxy'] = 'http://127.0.0.1:7890'
os.environ['https_proxy'] = 'http://127.0.0.1:7890'
P = {'http':'http://127.0.0.1:7890','https':'http://127.0.0.1:7890'}

# Top cancer types by sample count
PROJECTS = [
    # Top 6: already downloaded BRCA + 5 key cancer types
    "TCGA-LUAD",   # Lung adenocarcinoma (2nd largest)
    "TCGA-SKCM",   # Melanoma (matches our immunotherapy cohorts)
    "TCGA-COAD",   # Colorectal
    "TCGA-KIRC",   # Renal cell
    "TCGA-HNSC",   # Head and neck
]

BASE = Path("data/tcga")
BASE.mkdir(parents=True, exist_ok=True)

def gdc_query(endpoint, query_json, timeout=60):
    r = requests.post(f"https://api.gdc.cancer.gov/{endpoint}", json=query_json, proxies=P, verify=False, timeout=timeout)
    return r.json() if r.ok else None

def download_file(file_id, dest):
    if dest.exists(): return True
    for attempt in range(5):
        try:
            r = requests.get(f"https://api.gdc.cancer.gov/data/{file_id}", proxies=P, verify=False, timeout=300)
            if r.status_code == 200:
                with open(dest, 'wb') as f: f.write(r.content)
                return True
            time.sleep(5 * (attempt+1))
        except: time.sleep(10)
    return False

# 1. Download survival data
print("=== Downloading TCGA Survival Data ===")
surv_url = "https://api.gdc.cancer.gov/clinical_tar?size=11428"
r = requests.get(surv_url, proxies=P, verify=False, timeout=120)
if r.status_code == 200:
    surv_dest = BASE / "tcga_survival.tar.gz"
    with open(surv_dest, 'wb') as f: f.write(r.content)
    print(f"  Survival: {len(r.content)/1024/1024:.1f}MB")

# 2. Download expression for each project
total_mb = 0
for proj in PROJECTS:
    print(f"\n=== {proj} ===")
    proj_dir = BASE / proj
    proj_dir.mkdir(exist_ok=True)

    # Query files
    q = {
        "filters": {
            "op": "and",
            "content": [
                {"op":"in","content":{"field":"cases.project.project_id","value":[proj]}},
                {"op":"in","content":{"field":"files.data_type","value":["Gene Expression Quantification"]}},
            ]
        },
        "fields": "file_name,file_size,file_id,cases.case_id",
        "size": 1000
    }
    data = gdc_query("files", q)
    if not data:
        print(f"  Query failed"); continue

    hits = data['data']['hits']
    print(f"  {len(hits)} files")

    # Create case->file mapping, download one representative per case
    case_files = {}
    for hit in hits:
        for case in hit.get('cases', []):
            cid = case['case_id']
            if cid not in case_files or hit['file_size'] > case_files[cid][1]:
                case_files[cid] = (hit['file_id'], hit['file_size'], hit['file_name'])

    print(f"  {len(case_files)} unique cases")
    for i, (cid, (fid, size, fname)) in enumerate(case_files.items()):
        dest = proj_dir / f"{cid}.tsv.gz"
        if i < 3:  # Only show first few
            print(f"    {cid}: {size/1024/1024:.1f}MB - downloading...")
        if download_file(fid, dest):
            total_mb += dest.stat().st_size / 1024 / 1024
        if i % 10 == 0 and i > 0:
            print(f"    {i}/{len(case_files)} done ({total_mb:.0f}MB total)")
        time.sleep(0.5)  # Rate limit

print(f"\n{'='*60}")
print(f"TCGA download complete: {total_mb:.0f}MB total")
print(f"Projects: {len(PROJECTS)}, saved to {BASE}")
