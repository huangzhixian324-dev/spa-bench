"""v33 supplementary figures (S1-S12).

Regenerates every supplementary figure from the v33 authoritative outputs.
Output directory: results/figures/v33/ (the archived v28/v29 PNGs in
results/figures/ are left untouched).

Sources:
  results/benchmark/v33/benchmark_v33.json          (S1, S8, S12)
  results/benchmark/v33/permutation_v33.json        (S12)
  results/benchmark/v33/supplementary_v33_recompute.json
      (S2, S3, S4, S5, S6, S7, S10)
  data/tcga/tcga_cox_all.json                       (S11)
  manuscript Table 7 (power analysis)               (S9)

Figure S12 is the E4 full-permutation-matrix heatmap (strategy doc P1-1).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
V33 = REPO / "results" / "benchmark" / "v33"
OUTD = REPO / "results" / "figures" / "v33"
OUTD.mkdir(parents=True, exist_ok=True)

BM = json.load(open(V33 / "benchmark_v33.json"))["cohorts"]
PM = json.load(open(V33 / "permutation_v33.json"))["cohorts"]
RC = json.load(open(V33 / "supplementary_v33_recompute.json"))
MV = json.load(open(V33 / "mi_var_overlap_v33.json"))
TCGA = json.load(open(REPO / "data" / "tcga" / "tcga_cox_all.json"))

COHORT_LABELS = {
    "Hugo_2016": "Hugo 2016\nRECIST (n=28)",
    "Nathanson_2017": "Lauss 2017 (ACT)\nRECIST (n=25)",
    "Gide_2019_cBio": "Gide 2019\nRECIST (n=73)",
    "Jung_2019_DCB": "Jung 2019\nDCB (n=27)",
    "Riaz_2017_cytolytic": "Riaz 2017\ncytolytic (n=43)",
    "Riaz_2017_RECIST_v33": "Riaz 2017\nRECIST (n=42)",
}
METHODS = [("IMPRES", "IMPRES"), ("GEP", "GEP"), ("TIDE", "TIDE"),
           ("PD_L1", "PD-L1"), ("ElasticNet", "EN (MI)"),
           ("ElasticNet_Var", "EN (Var)")]
COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10,
                     "figure.dpi": 300, "savefig.bbox": "tight", "savefig.dpi": 300})


def save(fig, name):
    p = OUTD / name
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p.relative_to(REPO))


# ---------------------------------------------------------------- S1
def figS1():
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    x = np.arange(len(METHODS))
    cy = [BM["Riaz_2017_cytolytic"][m]["auroc"] for m, _ in METHODS]
    re = [BM["Riaz_2017_RECIST_v33"][m]["auroc"] for m, _ in METHODS]
    cy_ci = [BM["Riaz_2017_cytolytic"][m]["auroc_ci"] for m, _ in METHODS]
    re_ci = [BM["Riaz_2017_RECIST_v33"][m]["auroc_ci"] for m, _ in METHODS]
    cy_e = [[v - lo for v, (lo, _hi) in zip(cy, cy_ci)],
            [hi - v for v, (_lo, hi) in zip(cy, cy_ci)]]
    re_e = [[v - lo for v, (lo, _hi) in zip(re, re_ci)],
            [hi - v for v, (_lo, hi) in zip(re, re_ci)]]
    ax.bar(x - 0.2, cy, 0.38, label="Cytolytic endpoint (n=43)",
           color="#C44E52", yerr=cy_e, capsize=3,
           error_kw={"lw": 0.9, "ecolor": "#333333"})
    ax.bar(x + 0.2, re, 0.38, label="RECIST endpoint (n=42)",
           color="#4C72B0", yerr=re_e, capsize=3,
           error_kw={"lw": 0.9, "ecolor": "#333333"})
    for xi, v in zip(x - 0.2, cy):
        ax.text(xi, v + 0.02, f"{v:.3f}", ha="center", fontsize=7)
    for xi, v in zip(x + 0.2, re):
        ax.text(xi, v + 0.02, f"{v:.3f}", ha="center", fontsize=7)
    ax.axhline(0.5, ls="--", c="gray", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for _, lbl in METHODS])
    ax.set_ylabel("AUROC (leakage-free nested CV)")
    ax.set_ylim(0, 1.08)
    ax.set_title("Figure 1. Response-definition collapse (Riaz 2017, "
                 "HGNC-mapped, v33; 95% bootstrap CIs)")
    ax.legend(loc="upper left", fontsize=8)
    save(fig, "figS1_response_definition_collapse.png")


# ---------------------------------------------------------------- S2
def figS2():
    d = MV
    var = np.array(d["var_scores_all"])
    genes = d["genes_all"]
    gidx = {g: i for i, g in enumerate(genes)}
    mi_pre = {int(g): s for g, s in d["mi_prescreen_scores"]}  # idx -> MI
    top_var = set(d["top500_var_genes"])
    top_mi = set(d["top500_mi_prescreen_genes"])
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    pre_idx = np.array(sorted(mi_pre))
    v = var[pre_idx]
    mi_vals = np.array([mi_pre[i] for i in pre_idx])
    in_var = np.array([genes[i] in top_var for i in pre_idx])
    in_mi = np.array([genes[i] in top_mi for i in pre_idx])
    for m, c, lbl in [(~in_var & ~in_mi, "#BBBBBB", "prescreen, neither"),
                      (in_var & ~in_mi, "#55A868", "variance-only top-500"),
                      (~in_var & in_mi, "#4C72B0", "MI-only top-500"),
                      (in_var & in_mi, "#C44E52", "overlap (both)")]:
        ax.scatter(mi_vals[m], v[m], s=5, c=c, label=lbl, alpha=0.65, lw=0)
    ax.set_yscale("log")
    ax.set_xlabel("MI with response (stage-2, within top-2000 variance "
                  "prescreen; full data n=28)")
    ax.set_ylabel("Gene variance (log scale)")
    ax.set_title("Figure S2. MI vs. variance feature selection "
                 "(Hugo 2016, two-stage convention)")
    ax.legend(fontsize=7, loc="upper left")
    save(fig, "figS2_mi_vs_var.png")


# ---------------------------------------------------------------- S3
def figS3():
    from matplotlib.patches import Circle
    ov = MV["overlap_prescreen"]["count"]
    only = 500 - ov
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.add_patch(Circle((0.38, 0.5), 0.30, color="#4C72B0", alpha=0.45))
    ax.add_patch(Circle((0.62, 0.5), 0.30, color="#DD8452", alpha=0.45))
    ax.text(0.22, 0.5, f"MI only\n{only}", ha="center", va="center")
    ax.text(0.78, 0.5, f"Variance only\n{only}", ha="center", va="center")
    ax.text(0.5, 0.5, f"{ov}\n({100.0 * ov / 500:.0f}%)", ha="center",
            va="center", fontweight="bold")
    ax.text(0.5, 0.06, "Top-500 genes, Hugo 2016 (n=28): stage-2 MI within "
            "top-2000 variance prescreen (v33 pipeline convention)",
            ha="center", fontsize=8)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_title("Figure S3. Gene-set overlap: MI- vs. variance-selected")
    save(fig, "figS3_gene_overlap.png")


# ---------------------------------------------------------------- S4
def figS4():
    """True-event KM curves (iAtlas OS_STATUS), rebuilt 2026-09-07.

    Replaces the archived 0 < OS <= 400 d approximation figure. Curves
    are rebuilt from results/benchmark/v33/hugo_survival_true.json
    (per-patient iAtlas matches) plus the canonical IMPRES score
    recomputed from the local Hugo matrix; the resulting log-rank
    p-values are printed for verification against the archived JSON
    (IMPRES 0.554, RECIST 0.0003).
    """
    import re as _re
    import scanpy as sc
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    src = (REPO / "workflows" / "methods" / "base_method.py").read_text(
        encoding="utf-8")
    block = _re.search(r"PAIRS = \[(.*?)\]", src, _re.S).group(1)
    pairs = _re.findall(r'\("([A-Za-z0-9]+)",\s*"([A-Za-z0-9]+)"\)', block)

    d = json.load(open(V33 / "hugo_survival_true.json"))
    ad = sc.read_h5ad(REPO / "data" / "cohorts" / "Hugo_2016" /
                      "processed" / "Hugo_2016_processed.h5ad")
    X = ad.X.toarray() if hasattr(ad.X, "toarray") else np.asarray(ad.X)
    vidx = {g: i for i, g in enumerate(ad.var_names)}
    rows, resp, T, E = [], [], [], []
    for mt in d["matches"]:
        om = mt.get("os_months")
        if not om or (isinstance(om, float) and not np.isfinite(om)):
            continue
        rows.append(mt["local_row"])
        resp.append(mt["local_response"])
        T.append(om * 30.436875)
        E.append(1 if str(mt["os_status"]).startswith("1") else 0)
    scores = []
    for r in rows:
        hits = [1.0 if X[r, vidx[a]] > X[r, vidx[b]] else 0.0
                for a, b in pairs if a in vidx and b in vidx]
        scores.append(float(np.mean(hits)))
    scores = np.array(scores)
    T = np.array(T, dtype=float)
    E = np.array(E, dtype=float)
    resp = np.array(resp)
    hi = scores > np.median(scores)
    r1 = resp == 1
    p_imp = logrank_test(T[hi], T[~hi], E[hi], E[~hi]).p_value
    p_rec = logrank_test(T[r1], T[~r1], E[r1], E[~r1]).p_value
    print(f"[figS4] true-event log-rank: IMPRES p={p_imp:.4f} "
          f"(archived 0.554), RECIST p={p_rec:.6f} (archived 0.000326)")

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2))
    kmf = KaplanMeierFitter()
    ax = axes[0]
    for mask, col, lbl in ((hi, "#C44E52", f"IMPRES high (n={hi.sum()})"),
                           (~hi, "#4C72B0", f"IMPRES low (n={(~hi).sum()})")):
        kmf.fit(T[mask], E[mask], label=lbl)
        kmf.plot_survival_function(ax=ax, ci_show=False, color=col)
    ax.set_title(f"IMPRES median split (log-rank p = {p_imp:.3f})")
    ax.set_xlabel("Days"); ax.set_ylabel("Overall survival")
    ax.legend(fontsize=7)
    ax = axes[1]
    for mask, col, lbl in ((r1, "#55A868", "RECIST responder"),
                           (~r1, "#999999", "non-responder")):
        kmf.fit(T[mask], E[mask], label=lbl)
        kmf.plot_survival_function(ax=ax, ci_show=False, color=col)
    ax.set_title(f"RECIST v1.1 (log-rank p = {p_rec:.4f})")
    ax.set_xlabel("Days"); ax.set_ylabel("Overall survival")
    ax.legend(fontsize=7)
    fig.suptitle(f"Figure S4. Hugo 2016 overall survival (true iAtlas "
                 f"OS_STATUS events; n = {len(T)}, {int(E.sum())} events)",
                 y=1.02)
    save(fig, "figS4_km_survival.png")


# ---------------------------------------------------------------- S5
def figS5():
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.0))
    prep = RC["preprocessing_hugo"]
    labels, vals = [], []
    for k in ("fpkm_log2", "quantile"):
        if prep.get(k):
            labels.append(k.replace("_", " "))
            vals.append(prep[k]["auroc"])
    ax = axes[0]
    ax.bar(labels, vals, color=["#4C72B0", "#DD8452"], width=0.5)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=8)
    ax.axhline(0.5, ls="--", c="gray", lw=0.8)
    ax.set_ylim(0, 0.8); ax.set_ylabel("AUROC")
    ax.set_title("Preprocessing (Hugo, EN-MI)")
    cv = RC["cv_sensitivity_hugo"]
    keys = ["k=3", "k=5", "k=10", "LOO"]
    vals = [cv[k]["auroc"] for k in keys]
    ax = axes[1]
    ax.bar(keys, vals, color="#8172B3", width=0.55)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=8)
    ax.axhline(0.5, ls="--", c="gray", lw=0.8)
    ax.set_ylim(0, 0.8); ax.set_ylabel("AUROC")
    ax.set_title("Cross-validation strategy (Hugo, EN-MI)")
    fig.suptitle("Figure S5. Sensitivity analyses (v33 protocol; "
                 "TPM row not derivable from archived FPKM matrix)",
                 y=1.04)
    save(fig, "figS5_sensitivity.png")


# ---------------------------------------------------------------- S6
def figS6():
    rows = RC["learning_curve_hugo"]
    ns = [r["n"] for r in rows]
    mu = [r["auroc_mean"] for r in rows]
    sd = [r["auroc_sd_seeds"] for r in rows]
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    ax.errorbar(ns, mu, yerr=sd, marker="o", capsize=3, color="#4C72B0")
    ax.axhline(0.5, ls="--", c="gray", lw=0.8)
    ax.set_xlabel("Subsample size n (of 28)")
    ax.set_ylabel("AUROC (mean ± SD, 10 seeds)")
    ax.set_title("Figure S6. Learning curve (Hugo, ElasticNet-Var, "
                 "v33 protocol)")
    save(fig, "figS6_learning_curves.png")


# ---------------------------------------------------------------- S7
def figS7():
    d = RC["seed_stability_hugo"]
    dlt = np.array(d["mi_aurocs"]) - np.array(d["var_aurocs"])
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    ax.bar(np.arange(1, 11), dlt, color="#4C72B0", width=0.6)
    ax.axhline(d["delta_mean"], color="#C44E52",
               label=f"mean Δ = {d['delta_mean']:+.3f}")
    ax.axhspan(d["delta_mean"] - d["delta_sd"],
               d["delta_mean"] + d["delta_sd"], color="#C44E52", alpha=0.15,
               label=f"±1 SD = {d['delta_sd']:.3f}")
    ax.axhline(0, c="k", lw=0.8)
    ax.set_xlabel("Random seed")
    ax.set_ylabel("ΔAUROC (MI − Var)")
    ax.set_title("Figure S7. Seed stability (Hugo, 10 stratified seeds)")
    ax.legend(fontsize=7)
    save(fig, "figS7_seed_stability.png")


# ---------------------------------------------------------------- S8
def figS8():
    cohorts = list(COHORT_LABELS)
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 5.6))
    for ax, cname, col in zip(axes.ravel(), cohorts, COLORS):
        ys = np.arange(len(METHODS))[::-1]
        for y, (m, lbl), c in zip(ys, METHODS, COLORS):
            r = BM[cname][m]
            lo, hi = r["auroc_ci"]
            ax.plot([lo, hi], [y, y], color=c, lw=1.6)
            ax.plot(r["auroc"], y, "o", color=c, ms=5)
        ax.axvline(0.5, ls="--", c="gray", lw=0.8)
        ax.set_yticks(ys)
        ax.set_yticklabels([lbl for _, lbl in METHODS], fontsize=8)
        ax.set_xlim(0.0, 1.0)
        ax.set_title(COHORT_LABELS[cname].replace("\n", " "), fontsize=9)
        ax.set_xlabel("AUROC (95% bootstrap CI)", fontsize=8)
    fig.suptitle("Figure S8. Benchmark forest plot — 6 cohort-endpoints × "
                 "6 methods (v33)", y=1.0)
    fig.tight_layout()
    save(fig, "figS8_forest.png")


# ---------------------------------------------------------------- S9
def figS9():
    # Points reproduce manuscript Table 7 (asymptotic Mann-Whitney variance).
    rows = [(25, 0.40, 0.25), (28, 0.46, 0.22), (43, 0.21, 0.27),
            (73, 0.55, 0.15), (100, 0.50, 0.12), (200, 0.50, 0.08)]
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ns = [r[0] for r in rows]
    ds = [r[2] for r in rows]
    ax.plot(ns, ds, "o-", color="#4C72B0")
    for n, p, d in rows:
        ax.annotate(f"{int(p * 100)}%", (n, d), textcoords="offset points",
                    xytext=(4, 5), fontsize=7, color="#555555")
    ax.axhline(0.10, ls="--", c="gray", lw=0.8)
    ax.axhline(0.20, ls="--", c="gray", lw=0.8)
    ax.text(200, 0.10, "Δ = 0.10", fontsize=7, va="bottom", ha="right")
    ax.text(200, 0.20, "Δ = 0.20", fontsize=7, va="bottom", ha="right")
    ax.set_xlabel("Cohort size n")
    ax.set_ylabel("Min. detectable ΔAUROC at 80% power")
    ax.set_title("Figure S9. Power analysis (manuscript Table 8; % = responder "
                 "proportion)")
    save(fig, "figS9_power_curves.png")


# ---------------------------------------------------------------- S10
def figS10():
    cds = RC["cds_endpoints"]
    nulls = json.load(open(V33 / "cds_nulls_200_v33.json"))
    pts = []
    name_map = {"Hugo 2016 RECIST": ("Hugo_2016", "Hugo 2016"),
                "Nathanson 2017 RECIST": ("Nathanson_2017", "Lauss 2017 (ACT)"),
                "Gide 2019 RECIST": ("Gide_2019_cBio", "Gide 2019"),
                "Riaz 2017 RECIST (v33)": ("Riaz_2017_RECIST_v33",
                                           "Riaz 2017 RECIST"),
                "Riaz 2017 cytolytic": ("Riaz_2017_cytolytic",
                                        "Riaz 2017 cytolytic")}
    for label, (ck, disp) in name_map.items():
        v = cds.get(label)
        if not v:
            continue
        c = v.get("CDS", v.get("cds", v.get("total_score")))
        if c is None:
            continue
        ml = max(BM[ck]["ElasticNet"]["auroc"],
                 BM[ck]["ElasticNet_Var"]["auroc"])
        nl = nulls.get(label, {})
        pts.append((c, ml, disp, nl.get("mean"), nl.get("max")))
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.axvspan(0.70, 1.0, color="#C44E52", alpha=0.10)
    ax.axvline(0.70, ls="--", color="#C44E52", lw=1)
    ax.text(0.71, 0.33, "shipped HIGH boundary\n(advisory; Limitations)",
            fontsize=7.5, color="#C44E52")
    null_drawn = False
    for c, ml, lbl, nmean, nmax in pts:
        if nmean is not None:
            ax.plot([nmean, nmax], [ml - 0.035, ml - 0.035],
                    color="#555555", lw=5, alpha=0.35, zorder=1,
                    label=("cohort-matched null (mean\u2013max, 200 reps)"
                           if not null_drawn else None))
            null_drawn = True
        circ = "cytolytic" in lbl
        ax.scatter(c, ml, s=60, zorder=3,
                   color="#C44E52" if circ else "#4C72B0")
        ax.annotate(lbl, (c, ml), textcoords="offset points", xytext=(6, 4),
                    fontsize=8)
    ax.legend(loc="lower right", fontsize=7)
    ax.set_xlabel("CDS (v1.2.1, endpoint genes GZMA/PRF1)")
    ax.set_ylabel("Best trainable-ML AUROC (nested CV)")
    ax.set_xlim(0.3, 1.0); ax.set_ylim(0.3, 1.05)
    ax.set_title("Figure 2. CDS vs. observed ML performance (v33)")
    save(fig, "figS10_cds_vs_auroc.png")


# ---------------------------------------------------------------- S11
def figS11():
    sigs = sorted(TCGA["SKCM"], key=lambda s: TCGA["SKCM"][s]["HR"])
    meta = {"SKCM": (441, 212), "BRCA": (1082, 151), "LUAD": (510, 185),
            "COAD": (458, None)}
    cancers = [(c, f"{c} (n={TCGA[c][sigs[0]]['n']}, "
                   f"{TCGA[c][sigs[0]]['events']} events)")
               for c in ["SKCM", "BRCA", "LUAD", "COAD"] if c in TCGA]
    fig, axes = plt.subplots(1, len(cancers), figsize=(3.2 * len(cancers),
                                                        4.4),
                             sharey=True)
    if len(cancers) == 1:
        axes = [axes]
    ys = np.arange(len(sigs))[::-1]
    for ax, (ck, lbl) in zip(axes, cancers):
        for y, s in zip(ys, sigs):
            r = TCGA[ck][s]
            sig = r.get("p_bh", 1.0) < 0.05
            ax.plot(r["HR"], y, "o", ms=6,
                    color="#C44E52" if sig else "#BBBBBB",
                    fillstyle="full" if sig else "none")
        ax.axvline(1.0, ls="--", c="gray", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlabel("Hazard ratio (log scale)", fontsize=8)
        ax.set_title(lbl, fontsize=9)
    axes[0].set_yticks(ys)
    axes[0].set_yticklabels(sigs, fontsize=8)
    coad_note = ("; COAD re-acquired from GDC STAR counts" if "COAD" in TCGA
                 else "; COAD withdrawn")
    fig.suptitle(f"Figure S11. TCGA Cox per-SD HR, 12 immune signatures × "
                 f"{len(cancers)} cancer types ({coad_note}). Filled: BH "
                 f"q < 0.05 within cancer type.", y=1.06, fontsize=9)
    save(fig, "figS11_tcga_forest.png")


# ---------------------------------------------------------------- S12
def figS12():
    def bh(cells):
        # Standard BH adjusted p-values (cummin from largest rank).
        # 2026-09-07 fix: was cummax-from-smallest (overestimates q).
        lst = sorted(cells, key=lambda t: t[1])
        m, prev, q = len(lst), 1.0, {}
        for rank in range(m, 0, -1):
            k, p = lst[rank - 1]
            prev = min(prev, p * m / rank, 1.0)
            q[k] = prev
        return q

    cohorts = list(COHORT_LABELS)
    qv = {}
    for cname in cohorts:
        cells = [(m, PM[cname][m]["p"]) for m in PM[cname]]
        qv[cname] = bh(cells)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    grid = np.full((len(METHODS), len(cohorts)), np.nan)
    for i, (m, _) in enumerate(METHODS):
        for j, cname in enumerate(cohorts):
            cell = PM.get(cname, {}).get(m)
            if cell:
                grid[i, j] = -np.log10(max(cell["p"], 1e-6))
    im = ax.imshow(grid, cmap="Reds", vmin=0, vmax=5, aspect="auto")
    for i, (m, _) in enumerate(METHODS):
        for j, cname in enumerate(cohorts):
            cell = PM.get(cname, {}).get(m)
            if cell is None:
                note = ("excluded" if m == "PD_L1" else "n.c.")
                if m == "PD_L1" and cname:
                    note = "n.p."
                ax.text(j, i, note, ha="center", va="center", fontsize=7,
                        color="#888888")
                continue
            q = qv[cname][m]
            star = "*" if q < 0.05 else ""
            if m in ("ElasticNet", "ElasticNet_Var"):
                val = f"{BM[cname][m]['auroc']:.3f}/{cell['obs_auroc']:.3f}"
            else:
                val = f"{BM[cname][m]['auroc']:.3f}"
            ax.text(j, i - 0.16, f"{val}{star}",
                    ha="center", va="center", fontsize=6,
                    fontweight="bold" if star else "normal")
            ax.text(j, i + 0.22, f"p={cell['p']:.3g}", ha="center",
                    va="center", fontsize=6, color="#555555")
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels([COHORT_LABELS[c].replace("\n", " ") for c in cohorts],
                       fontsize=7, rotation=20, ha="right")
    ax.set_yticks(range(len(METHODS)))
    ax.set_yticklabels([lbl for _, lbl in METHODS], fontsize=8)
    plt.colorbar(im, ax=ax, label="−log10(permutation p)")
    ax.set_title("Figure 3. Complete permutation matrix (E4): AUROC + p per "
                 "cell; * BH q < 0.05 within cohort.\n"
                 "Fixed scorers: 50,000 prediction shuffles; trainable: "
                 "full-pipeline shuffles with frozen hyperparameters (per-cell "
                 "n); trainable cells show obs/perm-obs (p tests perm-obs)",
                 fontsize=8)
    save(fig, "figS12_permutation_heatmap.png")


def main():
    figS1(); figS2(); figS3(); figS4(); figS5()
    figS6(); figS7(); figS8(); figS9(); figS10()
    figS11(); figS12()
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
