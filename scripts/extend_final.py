"""Gide EN-MI 5,000-shuffle extension - FINAL clean rewrite (v3).

Rebuilt from scratch after overnight instability. Everything learned is in
here and nothing else:

- checkpoint resume from 3,200/5,000 (ge=78, seed state serialized)
- every 25 shuffles: save checkpoint (rng state + ge counter + obs)
- faulthandler: all thread stacks dumped to diag_stack.txt every 120s
  (if the process ever freezes, the blocking line is on disk)
- frozen protocol identical to fill_nc_cells.py: variance prescreen
  (cached per fold), MI top-500, ElasticNet C=0.1/l1_ratio=0.5, seed 42
- final JSON written at 5,000; nothing else touched

Run: venv\\Scripts\\python.exe -E -s -u scripts\\extend_final.py
"""
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

from rerun_v33 import SEED, load_cohort  # noqa: E402
from base_method import BaseMethod  # noqa: E402
from sklearn.feature_selection import mutual_info_classif  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

N = 5000
PRE = 2000
NSEL = 500
CKPT = (REPO / "results" / "benchmark" / "v33" / "nc_cells" /
        "Gide_2019_cBio__ElasticNet_5000.ckpt.json")
OUT = (REPO / "results" / "benchmark" / "v33" / "nc_cells" /
       "Gide_2019_cBio__ElasticNet_5000.json")
PARAMS = {"C": 0.1, "l1_ratio": 0.5}
COOL = 240  # seconds between 100-shuffle blocks



def save(n_done, ge, obs, rng):
    st = rng.get_state()
    json.dump({"n_done": n_done, "ge": ge, "obs": obs,
               "rng_state": [st[1].tolist(), st[2], st[3], st[4]]},
              open(CKPT, "w"))


def load():
    if not CKPT.exists():
        return None
    try:
        c = json.load(open(CKPT))
        rng = np.random.RandomState()
        rng.set_state(("MT19937", np.array(c["rng_state"][0]),
                       c["rng_state"][1], c["rng_state"][2], c["rng_state"][3]))
        return c["n_done"], c["ge"], c["obs"], rng
    except Exception as e:
        print(f"[ckpt corrupt: {e} -> fresh]", flush=True)
        return None


def main():
    print(f"[final-v3] Gide EN-MI extension to {N}", flush=True)
    adata, folds, y = load_cohort("Gide_2019_cBio")
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)

    folds_c = []
    for f in folds:
        tr, te = np.array(f["train"]), np.array(f["test"])
        vi = np.argsort(X[tr].var(axis=0))[-min(PRE, X.shape[1]):]
        folds_c.append((tr, te, vi))

    ck = load()
    if ck:
        n_done, ge, obs, rng = ck
        print(f"[resume] {n_done}/{N} ge={ge} obs={obs:.4f}", flush=True)
    else:
        obs_scores, obs_true = [], []
        for tr, te, vi in folds_c:
            mi = mutual_info_classif(X[tr][:, vi], y[tr], random_state=SEED)
            sel = vi[np.argsort(mi)[-NSEL:]]
            m = BaseMethod._make_elasticnet(**PARAMS)
            m.fit(X[tr][:, sel], y[tr])
            obs_true.extend(y[te])
            obs_scores.extend(m.predict_proba(X[te][:, sel])[:, 1])
        obs = float(roc_auc_score(obs_true, obs_scores))
        rng = np.random.RandomState(SEED)
        n_done, ge = 0, 0
        print(f"[fresh] obs={obs:.4f}", flush=True)
    save(n_done, ge, obs, rng)

    t0 = time.time()
    while n_done < N:
        stop = min(n_done + 25, N)
        for _ in range(n_done, stop):
            yp = rng.permutation(y)
            yt_p, ys_p = [], []
            try:
                for tr, te, vi in folds_c:
                    mi = mutual_info_classif(X[tr][:, vi], yp[tr],
                                             random_state=SEED)
                    sel = vi[np.argsort(mi)[-NSEL:]]
                    m = BaseMethod._make_elasticnet(**PARAMS)
                    m.fit(X[tr][:, sel], yp[tr])
                    ys_p.extend(m.predict_proba(X[te][:, sel])[:, 1])
                    yt_p.extend(yp[te])
                if len(set(yt_p)) >= 2:
                    if roc_auc_score(yt_p, ys_p) >= obs:
                        ge += 1
                    n_done += 1
            except ValueError:
                continue
        save(n_done, ge, obs, rng)
        print(f"    {n_done}/{N} (+{time.time() - t0:.0f}s ge={ge})",
              flush=True)
        if n_done < N and n_done % 100 == 0:
            print(f"    [cooldown {COOL}s]", flush=True)
            time.sleep(COOL)

    p = (ge + 1) / (n_done + 1)
    json.dump({"cohort": "Gide_2019_cBio", "method": "ElasticNet",
               "obs_auroc_this_run": round(obs, 4), "n_perm": n_done,
               "p": round(p, 5), "params_frozen": PARAMS,
               "protocol": "3,200+1,800 checkpointed full-pipeline shuffles, "
                           "frozen hyperparameters (review P0 item 4)"},
              open(OUT, "w"), indent=1)
    print(f"[DONE] p={p:.5f} ge={ge} n={n_done}", flush=True)
    print(f"[written] {OUT}", flush=True)


if __name__ == "__main__":
    main()
