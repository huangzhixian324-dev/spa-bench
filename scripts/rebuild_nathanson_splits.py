"""Rebuild the Lauss (GSE100797, 'Nathanson_2017') cohort from GEO FTP data.

Step 1 (this script): reconstruct the splits JSON and verify against the
run_registry fingerprint 12b146aea8c8d224. The split generator is not
archived, so we try the plausible deterministic constructions (stratified /
plain k-fold, shuffled or not, common seeds) against the recorded SHA-1.
"""
import gzip
import hashlib
import itertools
import json

import numpy as np
from sklearn.model_selection import StratifiedKFold, KFold

RAW = "data/cohorts/Nathanson_2017/raw"
TARGET = "12b146aea8c8d224"


def load_labels():
    titles, recists = [], []
    with gzip.open(f"{RAW}/GSE100797_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "recist" in line:
                recists = [t.strip('"').split(": ")[1]
                           for t in line.rstrip().split("\t")[1:]]
    y = [1 if r in ("CR", "PR") else 0 for r in recists]
    return titles, y


def dump_folds(folds, n, k):
    return json.dumps({"folds": folds, "n": n, "k": k})


def main():
    titles, y = load_labels()
    n = len(titles)
    print(f"n={n} samples, responders={sum(y)} ({sum(y)/n:.0%})")
    y_arr = np.array(y)

    for seed, shuffle, skf_cls in itertools.product(
            [42, 0, 1, 123, None], [True, False],
            [StratifiedKFold, KFold]):
        if skf_cls is KFold and shuffle is False:
            continue  # deterministic plain kfold equals no-shuffle variants
        kw = {"n_splits": 5}
        if shuffle:
            kw["random_state"] = seed if seed is not None else 42
            kw["shuffle"] = True
        else:
            kw["shuffle"] = False
        try:
            splitter = skf_cls(**kw)
            splits = splitter.split(np.zeros((n, 1)), y_arr)
            folds = []
            for i, (tr, te) in enumerate(splits):
                folds.append({"fold": i, "train": sorted(map(int, tr)),
                              "test": sorted(map(int, te))})
        except Exception:
            continue
        blob = dump_folds(folds, n, 5)
        h = hashlib.sha1(blob.encode()).hexdigest()[:16]
        tag = f"{skf_cls.__name__} seed={seed} shuffle={shuffle}"
        if h == TARGET:
            print(f"  *** MATCH  {tag}  sha1[:16]={h}")
            with open("data/cohorts/Nathanson_2017/processed/"
                      "Nathanson_2017_splits.json", "w") as fh:
                fh.write(blob)
            return 0
        print(f"  no    {tag}  {h}")
    print("no combination matched; need the original generator")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
