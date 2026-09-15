#!/usr/bin/env python3
"""
Robust GEO cohort downloader via HTTPS proxy with retry + SSL workaround.
"""
import os, re, time, requests, urllib3
from pathlib import Path
urllib3.disable_warnings()

os.environ.update({
    "http_proxy": "http://127.0.0.1:7890",
    "https_proxy": "http://127.0.0.1:7890",
})
PROXIES = {"http": "http://127.0.0.1:7890", "https": "http://127.0.0.1:7890"}

COHORTS = {
    "Hugo_2016": "GSE78220",
    "Gide_2019": "GSE115978",
    "Riaz_2017": "GSE91061",
    "Auslander_2018": "GSE119144",
    "Kim_2018": "GSE136961",
}


def robust_get(url, max_retries=5):
    """GET request with retry and backoff."""
    for i in range(max_retries):
        try:
            r = requests.get(url, proxies=PROXIES, verify=False, timeout=60)
            if r.status_code == 200:
                return r
            print(f"    HTTP {r.status_code}, retry {i+1}/{max_retries}")
        except Exception as e:
            print(f"    {type(e).__name__}, retry {i+1}/{max_retries}")
            time.sleep(2 * (i + 1))
    return None


def list_files(geo_id):
    """List supplementary files for a GEO series."""
    prefix = f"GSE{geo_id[3:5]}nnn"
    for pattern in [prefix, f"GSE{geo_id[:5]}nn"]:
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{pattern}/{geo_id}/suppl/"
        r = robust_get(url)
        if r and r.status_code == 200:
            files = re.findall(r'href=\"([^\"]+)\"', r.text)
            return [f for f in files if f not in ["/", "..", "/geo/"] and not f.startswith("http")]
    return []


def main():
    total = 0
    for name, geo_id in COHORTS.items():
        print(f"\n{'='*50}\n{name} ({geo_id})\n{'='*50}")
        out_dir = Path(f"data/cohorts/{name}/raw")
        out_dir.mkdir(parents=True, exist_ok=True)

        # 1. Download SOFT metadata
        prefix = f"GSE{geo_id[3:5]}nnn"
        soft_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/soft/{geo_id}_family.soft.gz"
        soft_dest = out_dir / f"{geo_id}_family.soft.gz"
        if not soft_dest.exists():
            r = robust_get(soft_url)
            if r:
                with open(soft_dest, "wb") as f:
                    f.write(r.content)
                print(f"  SOFT: {len(r.content)/1024:.0f} KB")

        # 2. List and download supplementary files
        files = list_files(geo_id)
        if not files:
            print(f"  No supplementary files found")
            continue

        for fname in files:
            dest = out_dir / fname
            if dest.exists():
                size = dest.stat().st_size / 1024 / 1024
                print(f"  {fname}: {size:.1f} MB (cached)")
                total += size
                continue

            file_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/suppl/{fname}"
            r = robust_get(file_url)
            if r:
                with open(dest, "wb") as f:
                    f.write(r.content)
                size = len(r.content) / 1024 / 1024
                print(f"  {fname}: {size:.1f} MB")
                total += size
            time.sleep(1)  # rate limit

    print(f"\n{'='*50}")
    print(f"Total downloaded: {total:.1f} MB across {len(COHORTS)} cohorts")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
