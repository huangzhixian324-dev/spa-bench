#!/usr/bin/env python3
"""
SPA-Bench Evaluation Metrics Module

Standardized evaluation metrics for immunotherapy response prediction.
All metrics computed in a unified way across all methods and cohorts.

Metrics:
  - Classification: AUROC, AUPRC, F1, Sensitivity, Specificity, Brier score
  - Survival: C-index, time-dependent AUROC
  - Clinical utility: Decision Curve Analysis (DCA), NRI, IDI
"""

import numpy as np
from scipy import stats
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score,
    recall_score, brier_score_loss, roc_curve, auc
)


# ============================================================
# Classification Metrics
# ============================================================

def classification_metrics(y_true, y_pred_prob, threshold=None):
    """
    Compute all classification metrics.

    Args:
        y_true: binary labels (0=non-responder, 1=responder)
        y_pred_prob: predicted probabilities
        threshold: classification threshold (None = 0.5)

    Returns:
        dict with all metrics
    """
    if threshold is None:
        threshold = 0.5

    y_pred_binary = (y_pred_prob >= threshold).astype(int)

    metrics = {
        "auroc": roc_auc_score(y_true, y_pred_prob),
        "auprc": average_precision_score(y_true, y_pred_prob),
        "f1": f1_score(y_true, y_pred_binary, zero_division=0),
        "sensitivity": recall_score(y_true, y_pred_binary, zero_division=0),
        "specificity": recall_score(1 - y_true, 1 - y_pred_binary, zero_division=0),
        "brier_score": brier_score_loss(y_true, y_pred_prob),
    }

    # Bootstrap 95% CI for AUROC
    metrics["auroc_ci_low"], metrics["auroc_ci_high"] = bootstrap_auroc_ci(
        y_true, y_pred_prob, n_bootstrap=1000, alpha=0.05
    )

    return metrics


def bootstrap_auroc_ci(y_true, y_pred_prob, n_bootstrap=1000, alpha=0.05):
    """Bootstrap 95% CI for AUROC."""
    n = len(y_true)
    aucs = []
    rng = np.random.RandomState(42)
    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_pred_prob[idx]))
    lower = np.percentile(aucs, 100 * alpha / 2)
    upper = np.percentile(aucs, 100 * (1 - alpha / 2))
    return lower, upper


# ============================================================
# Survival Metrics
# ============================================================

def c_index(y_true_time, y_true_event, y_pred_risk):
    """Harrell's C-index for survival data."""
    from scipy.stats import concordance
    result = concordance(y_true_time, y_pred_risk, y_true_event)
    return result.c_index


def time_dependent_auc(y_true_time, y_true_event, y_pred_risk,
                        eval_time=None):
    """Time-dependent AUROC at specified time point."""
    if eval_time is None:
        eval_time = np.median(y_true_time[y_true_event == 1])

    # Define cases/controls at eval_time
    is_case = (y_true_time <= eval_time) & (y_true_event == 1)
    is_control = y_true_time > eval_time

    y_binary = is_case.astype(int)
    mask = is_case | is_control

    if mask.sum() < 10:
        return np.nan

    return roc_auc_score(y_binary[mask], y_pred_risk[mask])


# ============================================================
# Clinical Utility Metrics
# ============================================================

def decision_curve_analysis(y_true, y_pred_prob, thresholds=None):
    """
    Decision Curve Analysis (Vickers & Elkin, 2006).

    Net Benefit = (TP - FP * w) / N
    where w = threshold / (1 - threshold)

    Returns:
        dict with thresholds and net_benefit for model, treat-all, treat-none
    """
    if thresholds is None:
        thresholds = np.linspace(0.01, 0.99, 99)

    n = len(y_true)
    net_benefit = []
    net_benefit_all = []
    net_benefit_none = np.zeros(len(thresholds))

    for t in thresholds:
        w = t / (1 - t)
        y_pred_binary = (y_pred_prob >= t).astype(int)
        tp = ((y_pred_binary == 1) & (y_true == 1)).sum()
        fp = ((y_pred_binary == 1) & (y_true == 0)).sum()
        nb = (tp - fp * w) / n
        net_benefit.append(nb)

        # Treat-all: everyone gets treatment
        nb_all = (y_true.sum() - (n - y_true.sum()) * w) / n
        net_benefit_all.append(nb_all)

    return {
        "thresholds": thresholds,
        "net_benefit": np.array(net_benefit),
        "net_benefit_all": np.array(net_benefit_all),
        "net_benefit_none": net_benefit_none,
    }

def _compute_nri(y_true, y_pred_old, y_pred_new, thresholds):
    """Internal: compute NRI components without bootstrap (avoids recursion)."""
    n_events = y_true.sum()
    n_nonevents = (1 - y_true).sum()

    def categorize(prob, thresholds):
        cat = np.zeros(len(prob), dtype=int)
        for i, t in enumerate(thresholds):
            cat[prob >= t] = i + 1
        return cat

    cat_old = categorize(y_pred_old, thresholds)
    cat_new = categorize(y_pred_new, thresholds)

    events_mask = y_true == 1
    if n_events > 0:
        nri_events = (cat_new[events_mask] > cat_old[events_mask]).sum() / n_events
        nri_events -= (cat_new[events_mask] < cat_old[events_mask]).sum() / n_events
    else:
        nri_events = 0.0

    nonevents_mask = y_true == 0
    nri_nonevents = (cat_new[nonevents_mask] < cat_old[nonevents_mask]).sum() / n_nonevents
    nri_nonevents -= (cat_new[nonevents_mask] > cat_old[nonevents_mask]).sum() / n_nonevents

    nri_total = nri_events + nri_nonevents
    return nri_events, nri_nonevents, nri_total


def net_reclassification_improvement(y_true, y_pred_old, y_pred_new,
                                      thresholds=None):
    """
    Net Reclassification Improvement (Pencina et al., 2008).

    NRI = NRI_events + NRI_nonevents
    """
    if thresholds is None:
        thresholds = [0.3, 0.5, 0.7]  # Low/Medium/High risk

    nri_events, nri_nonevents, nri_total = _compute_nri(
        y_true, y_pred_old, y_pred_new, thresholds)

    # Bootstrap test (uses _compute_nri internally, avoids recursion)
    nri_pval = _bootstrap_nri_pval(y_true, y_pred_old, y_pred_new,
                                   thresholds, n_bootstrap=500)

    return {
        "nri_events": nri_events,
        "nri_nonevents": nri_nonevents,
        "nri_total": nri_total,
        "nri_pval": nri_pval,
    }


def integrated_discrimination_improvement(y_true, y_pred_old, y_pred_new):
    """
    Integrated Discrimination Improvement (IDI).

    IDI = (IS_old - IS_new) where IS = integrated sensitivity.
    """
    # Mean predicted risk for events and non-events
    events = y_true == 1
    nonevents = y_true == 0

    idi_events = y_pred_new[events].mean() - y_pred_old[events].mean()
    idi_nonevents = y_pred_old[nonevents].mean() - y_pred_new[nonevents].mean()
    idi_total = idi_events + idi_nonevents

    return {
        "idi_events": idi_events,
        "idi_nonevents": idi_nonevents,
        "idi_total": idi_total,
    }


def _bootstrap_nri_pval(y_true, y_pred_old, y_pred_new, thresholds,
                         n_bootstrap=500):
    """Bootstrap test for NRI significance (internal, uses _compute_nri)."""
    rng = np.random.RandomState(42)
    n = len(y_true)
    nri_null = []
    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        _, _, nri_total = _compute_nri(
            y_true[idx], y_pred_old[idx], y_pred_new[idx], thresholds)
        nri_null.append(nri_total)

    # One-sided bootstrap p-value: fraction of bootstrap replicates whose
    # NRI is <= 0. NOTE (v33): this is the bootstrap distribution of the
    # observed statistic (centered near the observed NRI), not a permutation
    # null; interpret as a bootstrap-t style heuristic, and prefer the
    # permutation-based p-value for confirmatory claims.
    pval = (np.array(nri_null) <= 0).mean()
    return pval


# ============================================================
# Statistical Tests
# ============================================================

def delong_test(y_true, y_pred_a, y_pred_b):
    """Compare two AUROCs via a paired bootstrap test.

    NOTE (v33): despite the historical name, this is NOT the analytic DeLong
    covariance method — it is a paired bootstrap approximation (Sun & Xu 2014
    would be the fast exact alternative). The manuscript must refer to this
    test as a "paired bootstrap test of ΔAUROC", not "DeLong test".

    Returns:
        dict with statistic, p_value
    """
    n = len(y_true)
    # Compute AUC and its variance using DeLong method
    auc_a = roc_auc_score(y_true, y_pred_a)
    auc_b = roc_auc_score(y_true, y_pred_b)

    # Simplified DeLong: compute the variance-covariance
    # For a full implementation, see: Sun & Xu (2014)
    # Here we use a bootstrap approximation

    rng = np.random.RandomState(42)
    diffs = []
    for _ in range(1000):
        idx = rng.choice(n, size=n, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            continue
        diffs.append(
            roc_auc_score(y_true[idx], y_pred_a[idx]) -
            roc_auc_score(y_true[idx], y_pred_b[idx])
        )
    diffs = np.array(diffs)
    z_stat = np.mean(diffs) / (np.std(diffs) + 1e-10)
    p_val = 2 * stats.norm.sf(abs(z_stat))

    return {
        "auroc_a": auc_a,
        "auroc_b": auc_b,
        "delta_auroc": auc_a - auc_b,
        "z_statistic": z_stat,
        "p_value": p_val,
    }


# ============================================================
# Aggregate Eval
# ============================================================

def evaluate_all(y_true, y_pred_prob, cohort_name="unknown", y_pred_old=None):
    """Run all evaluations and return a structured report."""
    results = {"cohort": cohort_name, "n": len(y_true)}

    # Classification
    results.update(classification_metrics(y_true, y_pred_prob))

    # Clinical utility (vs baseline prediction)
    if y_pred_old is not None:
        dca = decision_curve_analysis(y_true, y_pred_prob)
        nri = net_reclassification_improvement(y_true, y_pred_old, y_pred_prob)
        idi = integrated_discrimination_improvement(y_true, y_pred_old, y_pred_prob)
        results["dca"] = dca
        results["nri"] = nri
        results["idi"] = idi

    return results
