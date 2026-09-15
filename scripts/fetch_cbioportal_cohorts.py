"""Fetch external ICI cohorts from cBioPortal for OS backfill + E3 backtest.

Cohorts (all reachable from this machine, unlike GEO):
  mel_iatlas_hugo_ucla_2016      Hugo 2016 (OS_STATUS truth backfill)
  blca_iatlas_imvigor210_2017    IMvigor210 (E3: the manuscript calls this
                                 'not publicly accessible' - it is, via
                                 iAtlas; CDS can now be computed)
  mel_iatlas_liu_2019            Liu 2019 DFCI anti-PD-1 melanoma (E3)
  mel_iatlas_riaz_nivolumab_2017 Riaz 2017 harmonized (OS for reference)
  skcm_vanderbilt_mskcc_2015     Van Allen 2015 samples (E3)

Uses POST /api/studies/export (official study exporter, returns the
original study files as a tar stream). Falls back to per-cohort tar if
the batch export fails.
"""
import io
import json
import sys
import tarfile
from pathlib import Path

import requests

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "external"
STUDIES = ["mel_iatlas_hugo_ucla_2016", "blca_iatlas_imvigor210_2017",
           "mel_iatlas_liu_2019", "mel_iatlas_riaz_nivolumab_2017",
           "skcm_vanderbilt_mskcc_2015"]
BASE = "https://www.cbioportal.org/api"


def export_studies(study_ids):
    r = requests.post(f"{BASE}/studies/export",
                      json={"studyIds": study_ids},
                      headers={"Accept": "application/x-tar"},
                      timeout=900, stream=True)
    r.raise_for_status()
    raw = r.content
    print(f"export payload: {len(raw) / 1e6:.1f} MB")
    tf = tarfile.open(fileobj=io.BytesIO(raw))
    n = 0
    for member in tf.getmembers():
        for sid in study_ids:
            if member.name.startswith(sid + "/") or member.name.startswith(sid + "\\"):
                tf.extract(member, OUT)
                n += 1
    print(f"extracted {n} members under {OUT}")
    return n


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        export_studies(STUDIES)
    except Exception as e:
        print(f"batch export failed ({e}); trying per-study", flush=True)
        for sid in STUDIES:
            try:
                export_studies([sid])
            except Exception as e2:
                print(f"[FAIL] {sid}: {e2}", flush=True)
    summary = {}
    for sid in STUDIES:
        d = OUT / sid
        files = sorted(p.name for p in d.rglob("*") if p.is_file())
        summary[sid] = files
        print(sid, "->", files[:6])
    json.dump(summary, open(OUT / "external_fetch_summary.json", "w"),
              indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
