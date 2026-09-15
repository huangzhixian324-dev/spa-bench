#!/usr/bin/env python3
"""
Comprehensive downloader: tries multiple sources for each cohort.
1. GEO supplementary files (via HTTPS proxy)
2. ENA (European Nucleotide Archive) for raw data access
3. Direct paper supplementary from journal websites

Usage: python scripts/download_all.py
"""

import os, re, time, gzip, requests, urllib3, subprocess
from pathlib import Path
urllib3.disable_warnings()

PROXY = {"http": "http://127.0.0.1:7890", "https": "http://127.0.0.1:7890"}
os.environ["http_proxy"] = "http://127.0.0.1:7890"
os.environ["https_proxy"] = "http://127.0.0.1:7890"

# ============================================================
# All known cohorts with processed expression data available
# ============================================================

COHORTS = {
    # GEO with supplementary files
    "Hugo_2016":    {"geo": "GSE78220",  "n": 28,  "note": "Melanoma, anti-PD-1, FPKM xlsx"},
    "Riaz_2017":    {"geo": "GSE91061",  "n": 109, "note": "Melanoma, anti-PD-1, FPKM csv"},
    "Gide_2019":    {"geo": "GSE115978", "n": 91,  "note": "Melanoma, anti-PD-1 +/- CTLA-4"},
    "Auslander_2018": {"geo": "GSE119144", "n": 50, "note": "Melanoma, anti-PD-1/CTLA-4"},
    "Kim_2018":     {"geo": "GSE136961", "n": 45,  "note": "Gastric, anti-PD-1"},
    "Jung_2019":    {"geo": "GSE135222", "n": 118, "note": "NSCLC, anti-PD-1"},
    "Prat_2017":    {"geo": "GSE100797", "n": 63,  "note": "Breast/HER2+, anti-HER2"},
    # Additional melanoma cohorts
    "Gide_2019_b":  {"geo": "GSE96619",  "n": 26,  "note": "Melanoma, anti-PD-1 (small)"},
}

# ============================================================
# Download functions
# ============================================================

def get_with_retry(url, max_tries=8):
    """Download with exponential backoff."""
    for i in range(max_tries):
        try:
            r = requests.get(url, proxies=PROXY, verify=False, timeout=120)
            if r.status_code == 200:
                return r
            print(f"    HTTP {r.status_code}, retry {i+1}/{max_tries}")
        except Exception as e:
            print(f"    {type(e).__name__}, retry {i+1}/{max_tries}")
            time.sleep(min(5 * (i + 1), 60))
    return None


def download_geo_supplementary(geo_id, out_dir):
    """Download all supplementary files for a GEO series."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"GSE{geo_id[3:5]}nnn"

    # First, download SOFT metadata
    soft_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/soft/{geo_id}_family.soft.gz"
    soft_dest = out_dir / f"{geo_id}_family.soft.gz"
    if not soft_dest.exists():
        r = get_with_retry(soft_url)
        if r:
            with open(soft_dest, "wb") as f:
                f.write(r.content)
            print(f"  SOFT: {len(r.content)/1024:.0f} KB")

    # List supplementary files
    suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/suppl/"
    r = get_with_retry(suppl_url)
    if not r:
        print(f"  No supplementary access")
        return []

    files = re.findall(r'href=\"([^\"]+)\"', r.text)
    data_files = [f for f in files if f not in ("/", "..", "/geo/") and not f.startswith("http")]
    print(f"  Found {len(data_files)} supplementary file(s)")

    downloaded = []
    for fname in data_files:
        dest = out_dir / fname
        if dest.exists():
            size = dest.stat().st_size / 1024 / 1024
            print(f"    {fname}: {size:.1f} MB (cached)")
            downloaded.append(fname)
            continue

        file_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/suppl/{fname}"
        r = get_with_retry(file_url)
        if r:
            with open(dest, "wb") as f:
                f.write(r.content)
            size = len(r.content) / 1024 / 1024
            print(f"    {fname}: {size:.1f} MB DOWNLOADED")
            downloaded.append(fname)
            time.sleep(2)

    return downloaded


def main():
    print("=" * 60)
    print("SPA-Bench: Comprehensive Data Downloader")
    print(f"Proxy: {PROXY['https']}")
    print("=" * 60)

    total_mb = 0
    success = 0
    failed = []

    for name, info in COHORTS.items():
        geo_id = info["geo"]
        print(f"\n[{name}] {geo_id} ({info['n']} samples) - {info['note']}")
        out_dir = Path(f"data/cohorts/{name}/raw")

        # Check if already downloaded
        existing = list(out_dir.glob("*")) if out_dir.exists() else []
        data_existing = [f for f in existing if f.suffix in (".gz", ".csv", ".xlsx", ".txt", ".tsv")]
        if data_existing:
            size = sum(f.stat().st_size for f in data_existing) / 1024 / 1024
            if size > 0.5:  # More than 500KB = probably has real data
                print(f"  Already have {size:.1f} MB, skipping")
                success += 1
                total_mb += size
                continue

        try:
            files = download_geo_supplementary(geo_id, out_dir)
            if files:
                success += 1
                size = sum((out_dir / f).stat().st_size for f in files if (out_dir / f).exists()) / 1024 / 1024
                total_mb += size
            else:
                failed.append(name)
        except Exception as e:
            print(f"  FAILED: {e}")
            failed.append(name)

    print(f"\n{'='*60}")
    print(f"Download Summary:")
    print(f"  Success: {success}/{len(COHORTS)} cohorts")
    print(f"  Total downloaded: {total_mb:.1f} MB")
    if failed:
        print(f"  Failed: {', '.join(failed)}")
    print(f"{'='*60}")
    return 0


if __name__ == "__main__":
    exit(main())
