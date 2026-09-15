#!/usr/bin/env python3
"""
Download publicly accessible immunotherapy cohort data from GEO.

Automatically downloads and extracts count matrices + clinical data
for cohorts that don't require dbGaP/EGA authentication.

Usage:
    python scripts/download_cohorts.py --all         # Download all GEO cohorts
    python scripts/download_cohorts.py --cohort Gide_2019
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

# ============================================================
# GEO-Accessible Cohorts (no auth required)
# ============================================================

GEO_COHORTS = {
    "Gide_2019": {
        "geo_id": "GSE115978",
        "cancer": "Melanoma",
        "n_samples": 100,
        "has_clinical": True,
        "clinical_file": "GSE115978_clinical.csv",
    },
    "Hugo_2016": {
        "geo_id": "GSE78220",
        "cancer": "Melanoma",
        "n_samples": 28,
        "has_clinical": True,
    },
    "Riaz_2017": {
        "geo_id": "GSE91061",
        "cancer": "Melanoma",
        "n_samples": 50,
        "has_clinical": True,
    },
    "Auslander_2018": {
        "geo_id": "GSE119144",
        "cancer": "Melanoma",
        "n_samples": 50,
        "has_clinical": True,
    },
    "Kim_2018": {
        "geo_id": "GSE136961",
        "cancer": "Gastric Cancer",
        "n_samples": 45,
        "has_clinical": True,
    },
}

# ============================================================
# EGA Cohorts (require pyega3 + credentials)
# ============================================================

EGA_COHORTS = {
    "IMvigor210": "EGAS00001002556",
    "OAK": "EGAD00001003610",
    "POPLAR": "EGAD00001003610",
}

# ============================================================
# dbGaP Cohorts (require dbGaP access)
# ============================================================

DBGAP_COHORTS = {
    "CheckMate_038": "phs001041",
    "Braun_2021": "phs002252",
    "VanAllen_2015": "phs000452",
    "Snyder_2017": "phs001178",
}


def download_geo_cohort(cohort_name, geo_id, output_dir):
    """Download count matrix from GEO using pysradb or direct curl."""
    output_dir = Path(output_dir) / cohort_name / "raw"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{cohort_name}] Downloading {geo_id} from GEO...")

    # Method 1: Use pysradb to get download links
    try:
        import pysradb
        sra = pysradb.SRAweb()
        df = sra.sra_metadata(geo_id)
        # Get run accessions
        runs = df['run_accession'].tolist()
        print(f"  Found {len(runs)} runs for {geo_id}")

        # Use fasterq-dump or prefetch
        for run in runs:
            print(f"  Downloading {run}...")
            subprocess.run([
                "prefetch", run, "-O", str(output_dir)
            ], check=False)
            subprocess.run([
                "fasterq-dump", run, "-O", str(output_dir),
                "--split-files"
            ], check=False)

    except ImportError:
        print("  pysradb not installed. Trying direct GEO download...")
        # Fallback: use GEOparse
        try:
            import GEOparse
            gse = GEOparse.get_GEO(geo=geo_id, destdir=str(output_dir),
                                   silent=True)
            print(f"  Downloaded {geo_id}: {len(gse.gsms)} samples")
        except ImportError:
            print(f"  [WARNING] Neither pysradb nor GEOparse available.")
            print(f"  [INFO] Please manually download {geo_id} from:")
            print(f"    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}")
            print(f"  [INFO] Save to: {output_dir}")

    print(f"[{cohort_name}] Download attempt complete.")


def main():
    parser = argparse.ArgumentParser(description="Download immunotherapy cohort data")
    parser.add_argument("--all", action="store_true",
                        help="Download all GEO-accessible cohorts")
    parser.add_argument("--cohort", help="Download a specific cohort")
    parser.add_argument("--output", default="data/cohorts",
                        help="Output directory")
    args = parser.parse_args()

    output_dir = Path(args.output)

    if args.all:
        print("=" * 60)
        print("Downloading all GEO-accessible cohorts")
        print(f"  Total: {len(GEO_COHORTS)} cohorts")
        print(f"  Output: {output_dir}")
        print("=" * 60)
        print()

        for name, info in GEO_COHORTS.items():
            download_geo_cohort(name, info["geo_id"], output_dir)

        print()
        print("=" * 60)
        print("GEO downloads complete!")
        print()
        print("Cohorts requiring additional access:")
        print("  EGA (pyega3 + credentials needed):")
        for name, acc in EGA_COHORTS.items():
            print(f"    - {name}: {acc}")
        print("  dbGaP (dbGaP access needed):")
        for name, acc in DBGAP_COHORTS.items():
            print(f"    - {name}: {acc}")
        print("=" * 60)

    elif args.cohort:
        if args.cohort in GEO_COHORTS:
            info = GEO_COHORTS[args.cohort]
            download_geo_cohort(args.cohort, info["geo_id"], output_dir)
        elif args.cohort in EGA_COHORTS:
            print(f"[{args.cohort}] Requires EGA authentication.")
            print(f"  pip install pyega3")
            print(f"  pyega3 fetch {EGA_COHORTS[args.cohort]} --output-dir {output_dir}/{args.cohort}/raw/")
        elif args.cohort in DBGAP_COHORTS:
            print(f"[{args.cohort}] Requires dbGaP access.")
        else:
            print(f"Unknown cohort: {args.cohort}")
            print(f"Available: {list(GEO_COHORTS.keys()) + list(EGA_COHORTS.keys()) + list(DBGAP_COHORTS.keys())}")

    else:
        parser.print_help()

    return 0


if __name__ == "__main__":
    exit(main())
