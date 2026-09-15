#!/usr/bin/env python3
"""
Download all GEO-accessible immunotherapy cohorts via proxy.

Usage:
    python scripts/download_all_cohorts.py
"""

import os
import re
import requests
from pathlib import Path

PROXY = {"http": "http://127.0.0.1:7890", "https": "http://127.0.0.1:7890"}
os.environ["http_proxy"] = "http://127.0.0.1:7890"
os.environ["https_proxy"] = "http://127.0.0.1:7890"

GEO_COHORTS = {
    "Hugo_2016": "GSE78220",
    "Gide_2019": "GSE115978",
    "Riaz_2017": "GSE91061",
    "Auslander_2018": "GSE119144",
    "Kim_2018": "GSE136961",
}


def list_suppl_files(geo_id):
    """List supplementary files for a GEO series."""
    prefix = f"GSE{geo_id[3:5]}nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/suppl/"
    try:
        resp = requests.get(url, proxies=PROXY, timeout=30)
        if resp.status_code != 200:
            return []
        files = re.findall(r'href=\"([^\"]+)\"', resp.text)
        return [f for f in files if f not in ["/", "..", "/geo/"] and
                not f.startswith("https://")]
    except Exception as e:
        print(f"  Error listing files: {e}")
        return []


def download_file(url, dest, desc=""):
    """Download a file with progress indicator."""
    try:
        resp = requests.get(url, proxies=PROXY, stream=True, timeout=120)
        if resp.status_code != 200:
            print(f"  {desc}: HTTP {resp.status_code}")
            return False
        total = int(resp.headers.get("content-length", 0))
        downloaded = 0
        with open(dest, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
        size_mb = downloaded / 1024 / 1024
        print(f"  {desc}: {size_mb:.1f} MB downloaded")
        return True
    except Exception as e:
        print(f"  {desc}: {e}")
        return False


def main():
    print("=" * 60)
    print("Downloading GEO Immunotherapy Cohorts via Proxy")
    print("=" * 60)

    for name, geo_id in GEO_COHORTS.items():
        print(f"\n[{name}] ({geo_id})")
        out_dir = Path(f"data/cohorts/{name}/raw")
        out_dir.mkdir(parents=True, exist_ok=True)

        # List available files
        files = list_suppl_files(geo_id)
        if not files:
            print(f"  No supplementary files found")
            continue

        print(f"  Found {len(files)} file(s): {files}")

        # Download each file
        prefix = f"GSE{geo_id[3:5]}nnn"
        base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/suppl"
        for f in files:
            url = f"{base_url}/{f}"
            dest = out_dir / f
            if dest.exists():
                print(f"  {f}: Already exists ({dest.stat().st_size/1024/1024:.1f} MB)")
                continue
            download_file(url, dest, f)

    print(f"\n{'=' * 60}")
    print("Download complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
