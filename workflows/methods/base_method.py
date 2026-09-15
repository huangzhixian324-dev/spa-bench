#!/usr/bin/env python3
"""
Base wrapper class for SPA-Bench method evaluation.

All methods (spatial, non-spatial, ML baselines) inherit from BaseMethod,
which enforces a standardized API for train/predict.
"""

import numpy as np
from abc import ABC, abstractmethod


class BaseMethod(ABC):
    """
    Abstract base class for all benchmarked methods.

    Every method must implement:
      - fit(adata, train_idx): Train on the specified samples
      - predict(adata, test_idx): Predict on the specified samples
    """

    def __init__(self, name, requires_spatial=False, requires_bulk=True):
        self.name = name
        self.requires_spatial = requires_spatial
        self.requires_bulk = requires_bulk

    @staticmethod
    def _to_dense(X):
        """v33: tolerate both sparse and dense AnnData.X throughout."""
        return X.toarray() if hasattr(X, "toarray") else np.asarray(X)

    @staticmethod
    def _make_elasticnet(C=1.0, l1_ratio=0.5, max_iter=10000):
        """Elastic-net logistic regression, version-adaptive.

        sklearn >= 1.8 deprecates penalty='elasticnet' and infers the penalty
        from l1_ratio instead. Earlier versions require penalty='elasticnet'
        explicitly, otherwise l1_ratio is silently ignored (which is exactly
        how the v28 code ended up fitting a plain L2 model). This helper
        guarantees elastic-net regularisation on both APIs.
        """
        from sklearn.linear_model import LogisticRegression
        import sklearn
        major, minor = (int(x) for x in sklearn.__version__.split(".")[:2])
        if (major, minor) >= (1, 8):
            return LogisticRegression(solver="saga", C=C, l1_ratio=l1_ratio,
                                      max_iter=max_iter)
        return LogisticRegression(solver="saga", penalty="elasticnet", C=C,
                                  l1_ratio=l1_ratio, max_iter=max_iter)

    @abstractmethod
    def fit(self, adata, train_idx):
        """Train the method on training samples."""
        pass

    @abstractmethod
    def predict(self, adata, test_idx):
        """Predict probabilities for test samples. Returns (N,) array."""
        pass

    def fit_predict(self, adata, train_idx, test_idx):
        """Combined fit + predict, with timing."""
        import time
        t0 = time.time()
        self.fit(adata, train_idx)
        fit_time = time.time() - t0

        t0 = time.time()
        preds = self.predict(adata, test_idx)
        pred_time = time.time() - t0

        return preds, {"fit_time": fit_time, "predict_time": pred_time}


# ============================================================
# Non-spatial baseline wrappers
# ============================================================

class PD_L1_Wrapper(BaseMethod):
    """PD-L1 IHC CPS/TPS as predictor."""
    def __init__(self):
        super().__init__("PD-L1_IHC", requires_spatial=False)

    def fit(self, adata, train_idx):
        pass  # No training needed

    def predict(self, adata, test_idx):
        # Use PD-L1 expression or clinical CPS score if available
        if "PDL1_CPS" in adata.obs.columns:
            scores = adata.obs["PDL1_CPS"].values[test_idx]
        elif "CD274" in adata.var_names:
            scores = self._to_dense(adata[test_idx, "CD274"].X).flatten()
        else:
            scores = np.zeros(len(test_idx))
        # Normalize to [0, 1]
        if scores.max() > scores.min():
            scores = (scores - scores.min()) / (scores.max() - scores.min())
        return scores


class TMB_Wrapper(BaseMethod):
    """Tumor Mutational Burden as predictor."""
    def __init__(self):
        super().__init__("TMB", requires_spatial=False)

    def fit(self, adata, train_idx):
        pass

    def predict(self, adata, test_idx):
        if "TMB" in adata.obs.columns:
            tmb = adata.obs["TMB"].values[test_idx].astype(float)
        else:
            tmb = np.ones(len(test_idx)) * 10  # default 10 muts/Mb
        # log-transform
        tmb = np.log2(np.clip(tmb, 1, None))
        if tmb.max() > tmb.min():
            tmb = (tmb - tmb.min()) / (tmb.max() - tmb.min())
        else:
            tmb = np.full(len(test_idx), 0.5)
        return tmb


class GEP_Wrapper(BaseMethod):
    """T-effector Gene Expression Profile signature."""
    def __init__(self):
        super().__init__("GEP", requires_spatial=False)
        self.genes = ["CD8A", "GZMA", "GZMB", "IFNG", "CXCL9",
                       "CXCL10", "PRF1", "TBX21"]

    def fit(self, adata, train_idx):
        pass

    def predict(self, adata, test_idx):
        available = [g for g in self.genes if g in adata.var_names]
        if not available:
            return np.zeros(len(test_idx))
        scores = self._to_dense(adata[test_idx, available].X).mean(axis=1)
        if scores.max() > scores.min():
            scores = (scores - scores.min()) / (scores.max() - scores.min())
        return scores

class TIDE_Wrapper(BaseMethod):
    """TIDE score (Jiang et al., Nature Medicine 2018) using official tidepy v1.3.9."""
    def __init__(self):
        super().__init__("TIDE", requires_spatial=False)
        self._tide_scores = None

    def fit(self, adata, train_idx):
        # TIDE is pre-trained; compute once on all data to avoid repeated calls.
        if self._tide_scores is None:
            try:
                from tidepy.pred import TIDE as tidepy_predict
                import pandas as pd
                X = adata.X
                X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
                # Reconstruct raw FPKM from log2(x+1) AnnData
                raw_X = np.exp2(X) - 1
                raw_X = np.clip(raw_X, 0, None)
                genes = adata.var_names.tolist()
                expr_df = pd.DataFrame(
                    raw_X.T, index=genes,
                    columns=[f"S{i}" for i in range(raw_X.shape[0])])
                result = tidepy_predict(expr_df, cancer="Melanoma", ignore_norm=True)
                self._tide_scores = -result["TIDE"].values  # negate
            except Exception as e:
                print(f"[WARN] Official TIDE unavailable: {e}")  # fallback below
                self._tide_scores = None

    def predict(self, adata, test_idx):
        if self._tide_scores is not None:
            scores = self._tide_scores[test_idx]
            s_min, s_max = scores.min(), scores.max()
            if s_max > s_min:
                scores = (scores - s_min) / (s_max - s_min)
            return scores
        # Fallback: simplified signature-based TIDE
        dys_genes = ["CXCL13","HAVCR2","PDCD1","TIGIT","LAG3","CTLA4"]
        exc_genes = ["VEGFA","VEGFB","TGFB1","CCL2","CXCL12","CXCL8"]
        avail_d = [g for g in dys_genes if g in adata.var_names]
        avail_e = [g for g in exc_genes if g in adata.var_names]
        d = self._to_dense(adata[test_idx, avail_d].X).mean(axis=1) if avail_d else np.zeros(len(test_idx))
        e = self._to_dense(adata[test_idx, avail_e].X).mean(axis=1) if avail_e else np.zeros(len(test_idx))
        z = np.clip(d - e - (d.mean() - e.mean()), -50, 50)
        return 1 / (1 + np.exp(z))



class IMPRES_Wrapper(BaseMethod):
    """IMPRES score (Auslander et al., Nature Medicine 2018) — canonical.

    v33 fix (Major-7): the previous implementation scored an arbitrary chain
    of 15 "pairs" (with one duplicated pair) via sign(e1-e2), which is NOT the
    published algorithm. The canonical IMPRES score counts how many of 15
    fixed pairwise logical relations are fulfilled in a sample:

        score = (# fulfilled relations) / 15

    The 15 relations were selected once by the original authors on
    neuroblastoma spontaneous-regression data and were never retrained on
    melanoma. The list below was decoded from the authors' FEATS.mat
    (github.com/noamaus/IMPRES-codes) and cross-checked against the 9 pairs
    named in the paper text (9/9 match). Following the paper, when only m of
    the 15 features are measured the score is linearly rescaled by 15/m.
    Higher scores predict response.
    """
    # (gene_A, gene_B) meaning: relation "expr(A) > expr(B)" is fulfilled
    PAIRS = [
        ("CD274",   "C10orf54"),   # PD-L1 > VISTA
        ("CD86",    "CD200"),      # CD86 > CD200
        ("CD40",    "CD274"),      # CD40 > PD-L1
        ("CD28",    "CD276"),      # CD28 > B7-H3
        ("CD40",    "CD28"),       # CD40 > CD28
        ("TNFRSF14", "CD86"),      # HVEM > CD86
        ("CD27",    "PDCD1"),      # CD27 > PD-1
        ("CD28",    "CD86"),       # CD28 > CD86
        ("CD40",    "CD80"),       # CD40 > CD80
        ("CD40",    "PDCD1"),      # CD40 > PD-1
        ("CD80",    "TNFSF9"),     # CD80 > 4-1BBL
        ("CD86",    "HAVCR2"),     # CD86 > TIM-3
        ("CD86",    "TNFSF4"),     # CD86 > OX40L
        ("CTLA4",   "TNFSF4"),     # CTLA-4 > OX40L
        ("PDCD1",   "TNFSF4"),     # PD-1 > OX40L
    ]

    def __init__(self):
        super().__init__("IMPRES", requires_spatial=False)

    def fit(self, adata, train_idx):
        pass  # IMPRES is trained once on neuroblastoma; no retraining here.

    def predict(self, adata, test_idx):
        X = self._to_dense(adata[test_idx].X)
        genes = list(adata.var_names)
        n_fulfilled = np.zeros(len(test_idx))
        n_measured = 0
        for g_a, g_b in self.PAIRS:
            if g_a in genes and g_b in genes:
                e_a = X[:, genes.index(g_a)]
                e_b = X[:, genes.index(g_b)]
                n_fulfilled += (e_a > e_b).astype(float)
                n_measured += 1
        if n_measured == 0:
            return np.full(len(test_idx), 0.5)
        # Paper's scaling: multiply back to the 0-15 scale when features
        # are missing, then divide by 15.
        scores = n_fulfilled * (15.0 / n_measured) / 15.0
        return np.clip(scores, 0, 1)


# ============================================================
# ML Baseline wrappers
# ============================================================

class ElasticNetWrapper(BaseMethod):
    """Elastic Net logistic regression on bulk RNA-seq (MI-based selection).

    v33 fix: penalty='elasticnet' was previously missing, so sklearn silently
    fell back to L2 and l1_ratio was ignored. Hyperparameters are now tuned
    via 3-fold inner cross-validation inside each training fold, matching the
    manuscript description.
    """
    N_SELECT = 500  # top-K genes by MI (aligned with ElasticNet_Var500 naming)

    def __init__(self):
        super().__init__("ElasticNet", requires_spatial=False)

    def _make_model(self):
        return self._make_elasticnet()

    def _select_features(self, X, y):
        """Outcome-aware feature selection (training fold only, no leakage)."""
        from sklearn.feature_selection import mutual_info_classif
        mi_scores = mutual_info_classif(X, y, random_state=42)
        n_select = min(self.N_SELECT, X.shape[1])
        return np.argsort(mi_scores)[-n_select:]

    def fit(self, adata, train_idx):
        from sklearn.model_selection import GridSearchCV
        X = self._to_dense(adata[train_idx].X)
        y = adata.obs["response"].values[train_idx].astype(int)

        self.top_genes = self._select_features(X, y)

        self.model = GridSearchCV(
            self._make_elasticnet(),
            param_grid={'C': [0.01, 0.1, 1.0], 'l1_ratio': [0.1, 0.5, 0.9]},
            cv=3, scoring='roc_auc')
        self.model.fit(X[:, self.top_genes], y)

    def predict(self, adata, test_idx):
        X = self._to_dense(adata[test_idx].X)
        return self.model.predict_proba(X[:, self.top_genes])[:, 1]


class ElasticNetVarWrapper(BaseMethod):
    """Elastic Net logistic regression with variance-based feature selection.

    v33: previously referenced in the manuscript and result JSONs
    (ElasticNet_Var500) but not implemented anywhere in the repository.
    Selects top-K genes by variance within the training fold, then fits the
    same elastic-net logistic regression with inner 3-fold CV.
    """
    N_SELECT = 500

    def __init__(self):
        super().__init__("ElasticNet_Var", requires_spatial=False)

    def fit(self, adata, train_idx):
        from sklearn.model_selection import GridSearchCV
        from sklearn.linear_model import LogisticRegression
        X = self._to_dense(adata[train_idx].X)
        y = adata.obs["response"].values[train_idx].astype(int)

        n_select = min(self.N_SELECT, X.shape[1])
        self.top_genes = np.argsort(X.var(axis=0))[-n_select:]

        self.model = GridSearchCV(
            self._make_elasticnet(),
            param_grid={'C': [0.01, 0.1, 1.0], 'l1_ratio': [0.1, 0.5, 0.9]},
            cv=3, scoring='roc_auc')
        self.model.fit(X[:, self.top_genes], y)

    def predict(self, adata, test_idx):
        X = self._to_dense(adata[test_idx].X)
        return self.model.predict_proba(X[:, self.top_genes])[:, 1]


class XGBoostWrapper(BaseMethod):
    """XGBoost on bulk RNA-seq."""
    def __init__(self):
        super().__init__("XGBoost", requires_spatial=False)

    def fit(self, adata, train_idx):
        from xgboost import XGBClassifier
        from sklearn.feature_selection import mutual_info_classif

        self.model = XGBClassifier(
            n_estimators=100, max_depth=3,
            learning_rate=0.1, subsample=0.8,
            random_state=42, eval_metric='logloss'
        )
        X = self._to_dense(adata[train_idx].X)
        y = adata.obs["response"].values[train_idx].astype(int)

        # Outcome-aware feature selection: mutual information
        mi_scores = mutual_info_classif(X, y, random_state=42)
        n_select = min(300, X.shape[1])
        self.top_genes = np.argsort(mi_scores)[-n_select:]

        self.model.fit(X[:, self.top_genes], y)

    def predict(self, adata, test_idx):
        X = self._to_dense(adata[test_idx].X)
        return self.model.predict_proba(X[:, self.top_genes])[:, 1]


# ============================================================
# Factory
# ============================================================

def get_method(method_name):
    """Factory: return method instance by name."""
    registry = {
        "PD_L1_IHC": PD_L1_Wrapper,
        "TMB": TMB_Wrapper,
        "GEP": GEP_Wrapper,
        "TIDE": TIDE_Wrapper,
        "IMPRES": IMPRES_Wrapper,
        "ElasticNet": ElasticNetWrapper,
        "ElasticNet_Var": ElasticNetVarWrapper,
        "XGBoost": XGBoostWrapper,
    }
    if method_name in registry:
        return registry[method_name]()
    else:
        raise ValueError(f"Unknown method: {method_name}. "
                         f"Available: {list(registry.keys())}")
