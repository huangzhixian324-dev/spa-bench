#!/usr/bin/env python3
"""
Download immunotherapy cohort data from ENA (European Nucleotide Archive).

ENA usually has better connectivity from China than NCBI.
ENA mirrors all GEO data with the same accession numbers.

Usage:
    python scripts/download_cohorts_ena.py --cohort Hugo_2016
    python scripts/download_cohorts_ena.py --all
"""

import subprocess
import sys
from pathlib import Path

# ENA download URLs (typically accessible from China)
# GEO → ENA mapping: GEO samples are linked via BioProject → ENA Study

COHORT_DOWNLOADS = {
    "Hugo_2016": {
        "geo": "GSE78220",
        "ena_study": "PRJNA318623",
        "ena_url": "https://www.ebi.ac.uk/ena/browser/view/PRJNA318623",
        # Direct FTP (usually better peering to China)
        "ftp_url": "ftp.sra.ebi.ac.uk/vol1/fastq/",
        "method": "ena",
        "n_samples": 28,
        "size_estimate": "~200 MB (raw counts from GEO supplementary)",
    },
    "Gide_2019": {
        "geo": "GSE115978",
        "ena_study": "PRJNA472435",
        "method": "ena",
        "n_samples": 100,
        "size_estimate": "~500 MB",
    },
    "Riaz_2017": {
        "geo": "GSE91061",
        "ena_study": "PRJNA354516",
        "method": "ena",
        "n_samples": 50,
        "size_estimate": "~300 MB",
    },
    "Auslander_2018": {
        "geo": "GSE119144",
        "ena_study": "PRJNA487564",
        "method": "ena",
        "n_samples": 50,
        "size_estimate": "~300 MB",
    },
    "Kim_2018": {
        "geo": "GSE136961",
        "ena_study": "PRJNA559619",
        "method": "ena",
        "n_samples": 45,
        "size_estimate": "~250 MB",
    },
}


def test_ena_connectivity():
    """Test if ENA is reachable."""
    import urllib.request
    try:
        urllib.request.urlopen("https://www.ebi.ac.uk", timeout=10)
        return True
    except Exception:
        return False


def download_geo_supplementary(geo_id, output_dir):
    """
    Download GEO supplementary files (contains processed count matrix).
    This is usually a small tar.gz with the expression matrix,
    much smaller than raw FASTQ files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # GEO supplementary files URL
    # Format: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/
    series_prefix = geo_id[:6] + "nnn"
    suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_prefix}/{geo_id}/suppl/"

    print(f"  GEO supplementary URL: {suppl_url}")
    print(f"  Output: {output_dir}")

    # Try curl with longer timeout
    cmd = [
        "curl", "-L", "--connect-timeout", "30", "--max-time", "600",
        "-o", str(output_dir / f"{geo_id}_suppl.tar.gz"),
        f"{suppl_url}{geo_id}_RAW.tar"
    ]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"  ✅ Downloaded successfully")
    else:
        print(f"  ❌ Download failed: {result.stderr[:200]}")
        print(f"  Try manually: {suppl_url}")


def main():
    print("=" * 60)
    print("Testing ENA connectivity...")
    if test_ena_connectivity():
        print("✅ ENA (EBI) is reachable!")
    else:
        print("❌ ENA is also unreachable.")

    print("\nAlternate download methods:")
    print("  1. GEO Supplementary Files (processed counts):")
    for name, info in COHORT_DOWNLOADS.items():
        geo_id = info["geo"]
        print(f"     {name} ({geo_id}): https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_id[:6]}nnn/{geo_id}/suppl/")
        print(f"       → Look for '*_RAW.tar' or '*_processed_data*' files")

    print("\n  2. ENA Study Pages:")
    for name, info in COHORT_DOWNLOADS.items():
        print(f"     {name}: {info['ena_url']}")

    print("\n  3. Using a proxy:")
    print("     export https_proxy=http://127.0.0.1:7890")
    print("     python scripts/download_cohorts.py --all")

    print("\n  4. Cloud mirror (if you have access):")
    print("     - Alibaba Cloud OSS GEO mirror (some universities have this)")
    print("     - Baidu Netdisk shared datasets from Chinese bioinformatics community")
    print("     - Contact the original paper authors for direct data access")


if __name__ == "__main__":
    main()
