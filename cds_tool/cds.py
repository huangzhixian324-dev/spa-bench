#!/usr/bin/env python3
"""
Circularity Detection Score (CDS) v1.0
=======================================
A tool for detecting circular dependencies between benchmark endpoints
and prediction features in transcriptomic studies.

Usage:
    # Command line
    cds --expr expression.csv --labels response.csv \
        --endpoint-genes GZMA,PRF1 --output report.html

    # Python API
    from cds import circularity_detection_score, CDS_REPORT_THRESHOLDS
    result = circularity_detection_score(X, y, ["GZMA","PRF1"], gene_list)
    print(f"CDS={result['CDS']:.3f} [{result['risk_level']}]")

Reference:
    SPATBench: A Standardized Framework Revealing Response Definition
    and Feature Selection as Critical Confounders. (2026)

Author: SPATBench Team
License: MIT
"""

import numpy as np
from sklearn.feature_selection import mutual_info_classif
from scipy.stats import spearmanr
import json, sys, argparse
from pathlib import Path

__version__ = "1.2.1"

# Risk thresholds
# v33 note: the shipped composite weights are C1 30% / C2 30% / C3 40% and
# these thresholds (HIGH > 0.70, MODERATE > 0.30). They were calibrated on a
# synthetic null of 30 shuffled-endpoint replicates AFTER inspecting the four
# real benchmark endpoints (disclosed in the manuscript); treat the risk
# labels as PROVISIONAL until external validation completes. Do not change
# these constants without recalibrating and re-versioning the tool.
#
# v1.2.0 note (2026-09-08): the categorical labels are screening heuristics
# with known calibration limits on real tumour RNA-seq — see CDS_CAVEATS.
# The validated reading is the CONTINUOUS value interpreted against a
# cohort-matched random-label null (pass null_values=... to enable the
# null-referenced output block). The primary circularity protection remains
# PROCEDURAL: declare the endpoint-defining genes and treat any
# gene-defined surrogate as a circular task by construction.
CDS_THRESHOLDS = {
    "HIGH": 0.70,
    "MODERATE": 0.30,
    "LOW": 0.00,
}

CDS_CAVEATS = [
    "Categorical labels (HIGH/MODERATE/LOW) are screening heuristics, "
    "not diagnoses; they have no specificity among structured biological "
    "endpoints (all 60 gene-program pseudo-endpoints in the SPATBench "
    "stress test scored HIGH).",
    "On real bulk tumour RNA-seq the C3 component typically saturates at "
    "1.0, bounding the composite to [0.40, 1.00]; the shipped LOW tier "
    "(<= 0.30) is therefore usually unreachable on such data.",
    "Under the shipped 0.70 HIGH boundary, 20-36% of 200 cohort-matched "
    "random-label nulls exceed it on real feature spaces (versus ~3% on "
    "the synthetic calibration null) - the shipped labels overstate "
    "specificity on real data.",
    "A HIGH reading on a CLINICAL endpoint can reflect genuine biology "
    "indistinguishable from circularity under CDS v1.x (e.g., Gide 2019: "
    "CDS 0.852 while GEP 0.830 is genuine signal). Always report the "
    "continuous value against a cohort-matched null.",
    "The primary circularity protection is procedural: declare the "
    "endpoint-defining genes and treat any gene-defined surrogate as a "
    "circular task by construction, whatever its score.",
]

def circularity_detection_score(X, y, endpoint_genes, gene_list,
                                 n_features=500, random_state=42,
                                 null_values=None, null_replicates=0,
                                 null_margin=0.1):
    """
    Compute Circularity Detection Score (CDS).

    CDS measures how much a gene-expression-derived response endpoint
    overlaps with features available for prediction.

    Components:
      C1 (30%): Mean MI between endpoint genes and response label
      C2 (30%): Mean |Spearman rho| between endpoint genes and response
      C3 (40%): MI skew ratio (top-feature MI / random-feature MI)

    Args:
        X: (n_samples, n_genes) expression matrix
        y: (n_samples,) binary response labels
        endpoint_genes: list of gene names/IDs used to define the endpoint
        gene_list: list of all gene names/IDs corresponding to X columns
        n_features: number of top-variance features to analyze (default 500)
        random_state: random seed for reproducibility

    Returns:
        dict with CDS, component scores, risk level, and interpretation.

    Interpretation:
        CDS > 0.70: HIGH circularity risk. Endpoint likely defined by genes in X.
                     ML methods will inflate AUROC. Use clinical endpoint instead.
        CDS 0.30-0.70: MODERATE. Some association; results should be interpreted cautiously.
        CDS < 0.30: LOW risk. Endpoint likely independent of feature space.

    Example:
        >>> import numpy as np
        >>> X = np.random.lognormal(4, 1, (100, 1000))
        >>> y = (X[:, 0] + X[:, 1] > np.median(X[:, 0] + X[:, 1])).astype(int)
        >>> genes = ["GENE_"+str(i) for i in range(1000)]
        >>> genes[0], genes[1] = "GZMA", "PRF1"
        >>> result = circularity_detection_score(X, y, ["GZMA","PRF1"], genes)
        >>> print(f"CDS={result['CDS']:.3f}, Risk={result['risk_level']}")
    """
    # C1: MI of endpoint genes with response
    ep_idx = [i for i, g in enumerate(gene_list) if g in endpoint_genes]
    if not ep_idx:
        return {
            "CDS": 0.0, "risk_level": "UNKNOWN",
            "error": f"None of {endpoint_genes} found in gene_list. "
                     f"Check gene identifier format.",
            "components": {"C1_MI": 0.0, "C2_rho": 0.0, "C3_skew": 0.0},
        }

    mi_vals = []
    for idx in ep_idx:
        mi = mutual_info_classif(X[:, idx:idx+1], y, random_state=random_state)
        mi_vals.append(mi[0])
    c1 = np.mean(mi_vals) if mi_vals else 0.0

    # C2: Spearman |rho| of endpoint genes with response
    rho_vals = []
    for idx in ep_idx:
        try:
            rho, _ = spearmanr(X[:, idx], y.astype(float))
            rho_vals.append(abs(rho))
        except:
            rho_vals.append(0.0)
    c2 = np.mean(rho_vals) if rho_vals else 0.0

    # C3: MI skew ratio (top-feature vs random-feature MI)
    n_top = min(n_features, X.shape[1])
    top_idx = np.argsort(X.var(axis=0))[-n_top:]
    mi_top = mutual_info_classif(X[:, top_idx], y, random_state=random_state)
    rng = np.random.RandomState(random_state)
    random_idx = rng.choice(X.shape[1], size=n_top, replace=False)
    mi_random = mutual_info_classif(X[:, random_idx], y, random_state=random_state)
    ratio = np.percentile(mi_top, 90) / max(np.median(mi_random), 0.0001)
    c3 = np.clip(ratio / 10.0, 0, 1)
    # v1.2.1: saturation marker. On real bulk RNA-seq the random-feature
    # median MI is ~0, so the denominator collapses and C3 clips to 1.0
    # regardless of endpoint circularity - a saturated C3 carries no
    # discrimination and the user should be told explicitly.
    c3_saturated = bool(ratio >= 10.0)
    c3_note = ("C3 is saturated (ratio/10 >= 1): the random-feature "
               "median MI is near zero, so the skew ratio hits the "
               "clipping ceiling. C3 carries no discrimination on this "
               "dataset; the composite is effectively "
               "0.4 + 0.3*C1n + 0.3*C2." if c3_saturated else "")

    # Composite CDS (weights: C1 30% / C2 30% / C3 40% — must match the
    # manuscript Methods; see v33 note at CDS_THRESHOLDS above)
    cds = 0.3 * np.clip(c1 / 0.05, 0, 1) + 0.3 * c2 + 0.4 * c3

    if cds > CDS_THRESHOLDS["HIGH"]:
        risk = "HIGH"
        guidance = (
            "CDS > 0.70 under the shipped screening heuristic. Two cases "
            "must be distinguished. (a) If the endpoint is a gene-defined "
            "surrogate, it is circular BY CONSTRUCTION (procedural rule): "
            "ML methods will learn the endpoint-defining genes, and the "
            "inflated AUROC is an artefact, not biology. (b) If the "
            "endpoint is a CLINICAL endpoint, the HIGH reading may reflect "
            "genuine biology that CDS v1.x cannot distinguish from "
            "circularity - report the continuous value against a "
            "cohort-matched random-label null and validate with an "
            "independent clinical endpoint; do NOT replace the clinical "
            "endpoint on this basis. See the caveats field for calibration "
            "limits of the shipped labels."
        )
    elif cds > CDS_THRESHOLDS["MODERATE"]:
        risk = "MODERATE"
        guidance = (
            "Some association between endpoint genes and prediction features exists "
            "(CDS 0.30-0.70). Results should be interpreted cautiously and stratified "
            "by CDS level. Consider sensitivity analyses with alternative endpoints. "
            "Interpret the continuous value against a cohort-matched random-label "
            "null rather than the categorical label alone."
        )
    else:
        risk = "LOW"
        guidance = (
            "CDS <= 0.30 under the shipped screening heuristic. Note that "
            "on real bulk tumour RNA-seq this tier is usually unreachable "
            "because C3 saturates at 1.0 (composite bounded to [0.40, "
            "1.00]); a LOW reading here should be double-checked against "
            "the C3 component value. Interpret against a cohort-matched "
            "random-label null."
        )

    result = {
        "CDS": round(float(cds), 3),
        "risk_level": risk,
        "guidance": guidance,
        "components": {
            "C1_MI_endpoint_genes": round(float(c1), 4),
            "C2_rho_endpoint_response": round(float(c2), 4),
            "C3_MI_skew_ratio": round(float(c3), 4),
            "C3_saturated": c3_saturated,
        },
        "endpoint_genes_found": [gene_list[i] for i in ep_idx],
        "version": __version__,
        "caveats": list(CDS_CAVEATS),
    }
    if c3_note:
        result["caveats"] = [c3_note] + result["caveats"]

    # ---- built-in null generator (v1.2.1): the validated reading is one
    # call away. Explicit null_values take precedence over generated nulls.
    if null_values is None and null_replicates > 0:
        rng_null = np.random.RandomState(random_state)
        null_values = []
        for _ in range(null_replicates):
            y_perm = rng_null.permutation(np.asarray(y))
            mi_t = mutual_info_classif(X[:, top_idx], y_perm,
                                       random_state=random_state)
            mi_r = mutual_info_classif(X[:, random_idx], y_perm,
                                       random_state=random_state)
            r_perm = (np.percentile(mi_t, 90) /
                      max(np.median(mi_r), 0.0001))
            c3_p = np.clip(r_perm / 10.0, 0, 1)
            mi_ep_p = []
            for idx in ep_idx:
                v = mutual_info_classif(X[:, idx:idx + 1], y_perm,
                                        random_state=random_state)
                mi_ep_p.append(v[0])
            c1_p = np.mean(mi_ep_p) if mi_ep_p else 0.0
            rho_p = []
            for idx in ep_idx:
                try:
                    rp, _ = spearmanr(X[:, idx], y_perm.astype(float))
                    rho_p.append(abs(rp))
                except Exception:
                    rho_p.append(0.0)
            c2_p = np.mean(rho_p) if rho_p else 0.0
            null_values.append(
                0.3 * np.clip(c1_p / 0.05, 0, 1) + 0.3 * c2_p + 0.4 * c3_p)

    # ---- null-referenced reading (validated claim of the SPATBench paper)
    if null_values is not None:
        nv = np.asarray(null_values, dtype=float)
        nv = nv[np.isfinite(nv)]
        null_mean = float(np.mean(nv))
        null_max = float(np.max(nv))
        if float(cds) > null_max + null_margin:
            reading = "far above null max"
        elif float(cds) > null_max:
            reading = "above null max"
        elif float(cds) > null_mean:
            reading = "above null mean"
        else:
            reading = "at or below null mean"
        result["null_referenced"] = {
            "n_null_replicates": int(nv.size),
            "null_mean": round(null_mean, 4),
            "null_max": round(null_max, 4),
            "reading": reading,
            "note": "Validated reading mode: interpret the continuous CDS "
                    "against the cohort-matched random-label null. A HIGH "
                    "categorical label with a clinical endpoint does not "
                    "establish circularity (see caveats).",
        }
        # prepend the null-referenced reading to guidance for CLI users
        result["guidance"] = (
            f"Null-referenced reading: CDS {float(cds):.3f} is '{reading}' "
            f"(cohort-matched null: mean {null_mean:.3f}, max {null_max:.3f}, "
            f"{nv.size} replicates). " + result["guidance"])
    return result


def batch_cds(expression_file, labels, endpoint_genes, gene_id_col=0, sep=None,
              null_values=None, null_replicates=0):
    """
    Convenience function: load data from files and compute CDS.

    Args:
        expression_file: path to CSV/TSV expression matrix (genes x samples)
        labels: path to CSV/TSV with binary response labels
        endpoint_genes: list of gene names
        gene_id_col: column index for gene IDs (default 0 = first column)

    Returns:
        CDS result dict
    """
    import pandas as pd
    if sep is None:
        sep = "\t" if expression_file.endswith(".tsv") else ","
    df = pd.read_csv(expression_file, sep=sep, index_col=gene_id_col)
    lbl = pd.read_csv(labels, sep=sep, index_col=0)
    X = df.T.values
    y = lbl.iloc[:, 0].values.astype(int)
    gene_list = df.index.tolist()
    return circularity_detection_score(X, y, endpoint_genes, gene_list,
                                       null_values=null_values,
                                       null_replicates=null_replicates)


def generate_html_report(result, output_path=None):
    """Generate a self-contained HTML report from CDS results."""
    risk_color = {"HIGH": "#e65100", "MODERATE": "#f57c00", "LOW": "#2e7d32", "UNKNOWN": "#999"}
    color = risk_color.get(result.get("risk_level", "UNKNOWN"), "#999")

    nr = result.get("null_referenced")
    nr_html = ""
    if nr:
        nr_html = f"""
<div class="card">
<h2>Null-Referenced Reading (validated)</h2>
<p>CDS <strong>{result['CDS']:.3f}</strong> is <strong>'{nr['reading']}'</strong> against the
cohort-matched random-label null (mean {nr['null_mean']:.3f}, max {nr['null_max']:.3f},
{nr['n_null_replicates']} replicates).</p>
<p style="color:#666;font-size:0.92em">{nr.get('note','')}</p>
</div>"""

    caveats_html = ""
    if result.get("caveats"):
        items = "".join(f"<li>{c}</li>" for c in result["caveats"])
        caveats_html = f"""
<div class="card">
<h2>Calibration Caveats (read before acting on the label)</h2>
<ul style="margin:8px 0 0 18px;color:#555;font-size:0.92em">{items}</ul>
</div>"""

    c3_sat = result["components"].get("C3_saturated")
    c3_badge = (' <span style="color:#c62828;font-weight:bold">'
                '[SATURATED - no discrimination]</span>' if c3_sat else "")

    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>CDS Report</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,sans-serif;max-width:700px;margin:40px auto;padding:20px;background:#fafafa}}
.card{{background:white;border-radius:8px;padding:24px;margin:16px 0;box-shadow:0 2px 8px rgba(0,0,0,0.08)}}
h1{{color:#1a237e}} h2{{color:#333;border-bottom:2px solid #e0e0e0;padding-bottom:8px}}
.big-score{{font-size:3em;font-weight:bold;text-align:center;color:{color};margin:20px 0}}
.risk-badge{{display:inline-block;padding:6px 16px;border-radius:20px;font-weight:bold;color:white;background:{color}}}
.guidance{{background:#fff3e0;border-left:4px solid {color};padding:16px;margin:16px 0;border-radius:4px}}
table{{width:100%;border-collapse:collapse;margin:12px 0}}
th,td{{padding:10px;text-align:left;border-bottom:1px solid #e0e0e0}}
th{{background:#f5f5f5}}
footer{{text-align:center;color:#999;font-size:0.85em;margin-top:30px}}
</style></head><body>
<h1>Circularity Detection Score (CDS) Report</h1>
<div class="card">
<h2>Result</h2>
<div class="big-score">{result['CDS']:.3f}</div>
<p style="text-align:center"><span class="risk-badge">{result['risk_level']} RISK</span></p>
<div class="guidance">{result.get('guidance', '')}</div>
</div>
{nr_html}
<div class="card">
<h2>Component Breakdown</h2>
<table>
<tr><th>Component</th><th>Value</th><th>Description</th></tr>
<tr><td>C1 (MI)</td><td>{result['components']['C1_MI_endpoint_genes']:.4f}</td><td>Mutual info: endpoint genes with response</td></tr>
<tr><td>C2 (rho)</td><td>{result['components']['C2_rho_endpoint_response']:.4f}</td><td>Spearman |rho|: endpoint genes with response</td></tr>
<tr><td>C3 (skew)</td><td>{result['components']['C3_MI_skew_ratio']:.4f}{c3_badge}</td><td>MI skew ratio: top features / random features</td></tr>
</table>
<p><strong>Endpoint genes found:</strong> {', '.join(result.get('endpoint_genes_found', ['none']))}</p>
</div>
{caveats_html}
<footer>CDS v{__version__} | SPATBench Project | MIT License</footer>
</body></html>"""

    if output_path:
        with open(output_path, "w") as f:
            f.write(html)
        return output_path
    return html


# ============================================================
# CLI
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Circularity Detection Score (CDS) ? detect circular dependencies in benchmark endpoints.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  cds --expr expression.tsv --labels response.csv \
      --endpoint-genes GZMA,PRF1

  # With HTML report
  cds --expr data.tsv --labels labels.csv \
      --endpoint-genes GZMA,PRF1 --html report.html

  # Demo with synthetic data
  cds --demo
        """
    )
    parser.add_argument("--expr", help="Expression matrix file (genes x samples, CSV/TSV)")
    parser.add_argument("--labels", help="Response labels file (CSV/TSV, first column = binary labels)")
    parser.add_argument("--endpoint-genes", help="Comma-separated endpoint-defining gene names")
    parser.add_argument("--null-replicates", type=int, default=0,
                        help="Generate this many cohort-matched random-label null replicates "
                             "and enable the null-referenced reading (recommended: 100-200)")
    parser.add_argument("--null-file", help="Optional JSON/CSV file with pre-computed null values "
                                            "(takes precedence over --null-replicates)")
    parser.add_argument("--html", help="Output HTML report path")
    parser.add_argument("--json", help="Output JSON results path")
    parser.add_argument("--demo", action="store_true", help="Run demo with synthetic data")
    parser.add_argument("--version", action="version", version=f"cds v{__version__}")

    args = parser.parse_args()

    null_values = None
    if args.null_file:
        nv = json.load(open(args.null_file))
        if isinstance(nv, dict) and "values" in nv:
            nv = nv["values"]
        null_values = nv

    if args.demo:
        print("CDS Demo: Synthetic data")
        print("=" * 50)
        np.random.seed(42)
        N, G = 200, 5000
        X = np.random.lognormal(4, 1, (N, G))
        genes = [f"GENE_{i}" for i in range(G)]
        genes[0], genes[1] = "GZMA", "PRF1"

        # Safe endpoint
        y_safe = np.random.binomial(1, 0.4, N)
        r_safe = circularity_detection_score(X, y_safe, ["GZMA","PRF1"], genes)
        print(f"\nSafe endpoint (random response):")
        print(f"  CDS = {r_safe['CDS']:.3f} [{r_safe['risk_level']}]")
        print(f"  C1={r_safe['components']['C1_MI_endpoint_genes']:.4f} "
              f"C2={r_safe['components']['C2_rho_endpoint_response']:.4f} "
              f"C3={r_safe['components']['C3_MI_skew_ratio']:.4f}")

        # Circular endpoint (with built-in null: the validated reading)
        y_circ = (X[:,0] + X[:,1] > np.median(X[:,0] + X[:,1])).astype(int)
        r_circ = circularity_detection_score(X, y_circ, ["GZMA","PRF1"], genes,
                                             null_replicates=50)
        print(f"\nCircular endpoint (GZMA+PRF1-defined, 50 built-in null replicates):")
        print(f"  CDS = {r_circ['CDS']:.3f} [{r_circ['risk_level']}]")
        print(f"  C1={r_circ['components']['C1_MI_endpoint_genes']:.4f} "
              f"C2={r_circ['components']['C2_rho_endpoint_response']:.4f} "
              f"C3={r_circ['components']['C3_MI_skew_ratio']:.4f} "
              f"saturated={r_circ['components']['C3_saturated']}")
        nr = r_circ.get('null_referenced')
        if nr:
            print(f"  Null-referenced: '{nr['reading']}' "
                  f"(null mean {nr['null_mean']:.3f}, max {nr['null_max']:.3f})")
        print(f"\n  Guidance: {r_circ['guidance'][:220]}...")

        if args.html:
            generate_html_report(r_circ, args.html)
            print(f"\nHTML report saved: {args.html}")
        return 0

    # Real data mode
    if not args.expr or not args.labels or not args.endpoint_genes:
        parser.error("--expr, --labels, and --endpoint-genes are required (or use --demo)")

    endpoint_genes = [g.strip() for g in args.endpoint_genes.split(",")]

    print(f"Loading: {args.expr}")
    result = batch_cds(args.expr, args.labels, endpoint_genes,
                       null_values=null_values,
                       null_replicates=args.null_replicates)

    print(f"\nCDS = {result['CDS']:.3f} [{result['risk_level']} RISK]")
    for k, v in result['components'].items():
        print(f"  {k}: {v:.4f}")
    print(f"\n{result.get('guidance', '')}")

    if args.html:
        generate_html_report(result, args.html)
        print(f"\nHTML report: {args.html}")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(result, f, indent=2)
        print(f"JSON: {args.json}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
