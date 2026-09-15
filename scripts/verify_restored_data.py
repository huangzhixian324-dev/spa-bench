"""Final equivalence verification: run the four fixed scorers (IMPRES/GEP/
TIDE/PD-L1) OOF AUROC on the RESTORED original h5ad+splits for Jung, Lauss
(Nathanson) and Gide, and compare against benchmark_v33.json.

If all 12 fingerprints match, the restored files ARE the original data.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
os.chdir(REPO)
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

from rerun_v33 import COHORTS, load_cohort  # noqa: E402
from base_method import get_method  # noqa: E402

BM = json.load(open(REPO / "results/benchmark/v33/benchmark_v33.json"))["cohorts"]
METHOD_MAP = [("IMPRES", "IMPRES"), ("GEP", "GEP"), ("TIDE", "TIDE"),
              ("PD_L1_IHC", "PD_L1")]
COHORT_KEYS = {"Jung_2019_DCB": "Jung_2019_DCB",
               "Nathanson_2017": "Nathanson_2017",
               "Gide_2019_cBio": "Gide_2019_cBio"}

all_ok = True
results = {}
for cohort, bm_key in COHORT_KEYS.items():
    print(f"=== {cohort} ===", flush=True)
    adata, folds, y = load_cohort(cohort)
    for meth, bm_meth in METHOD_MAP:
        m = get_method(meth)
        yt, yp = [], []
        for f in folds:
            preds, _ = m.fit_predict(adata, f["train"], f["test"])
            yt.extend(np.asarray(y)[f["test"]])
            yp.extend(preds)
        a = float(roc_auc_score(yt, yp))
        expect = BM[bm_key][bm_meth]["auroc"]
        ok = abs(a - expect) < 0.0005
        all_ok = all_ok and ok
        results[f"{cohort}/{bm_meth}"] = {
            "rebuilt": round(a, 4), "benchmark": round(expect, 4),
            "match": ok}
        print(f"  {meth:10s} rebuilt {a:.4f} | benchmark {expect:.4f} | "
              f"{'*** MATCH' if ok else 'DIFF ' + format(abs(a-expect), '.4f')}",
              flush=True)

print()
print("=" * 60)
print("FINAL:", "ALL 12 FINGERPRINTS MATCH - restored data IS the original data"
      if all_ok else "SOME FINGERPRINTS DIFFER - see above")
print("=" * 60)
json.dump(results, open(REPO / "results/benchmark/v33/"
                        "restore_fingerprint_check.json", "w"), indent=1)
print("saved restore_fingerprint_check.json")
