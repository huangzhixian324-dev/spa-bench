"""E2: CDS negative-control stress test (strategy doc P0-2).

Question: does CDS v1.1.0 distinguish gene-program-defined pseudo-endpoints
(biologically real but non-clinical) from genuinely circular endpoints and
from random labels?

Design:
  - 30 pseudo-endpoints = 3 families x 10 canonical gene programs
    (A: ECM/stroma without EMT-TFs; B: tumor purity/proliferation;
     C: lineage/metabolic/immune programs). No program is a clinical
     endpoint, and every program is REAL (canonical markers), unlike the
    purely random-label null.
  - Each pseudo-endpoint = median split of the program's mean
    log2(FPKM+1) expression (the same convention as the cytolytic
    surrogate), so every pseudo-endpoint IS gene-defined by construction.
  - Declared endpoint genes = the program's own genes (the intended CDS
    usage).
  - Cohorts: Gide 2019 (n = 73, largest) and Hugo 2016 (n = 28, small-n
    regime).
  - References computed identically on the same cohorts: the cytolytic
    surrogate (GZMA/PRF1 mean, median split - the true circular positive
    control) and 10 random-label nulls per cohort (seed 42).

Reading of the result (both interpretations are reported; neither is
spun): HIGH on pseudo-endpoints = correct detection of gene-defined
surrogates regardless of program identity (program-agnostic sensitivity);
agreement of pseudo-endpoint CDS with the random-label floor would instead
indicate the C1 low bar (MI >= 0.05 saturates C1) makes CDS flag any
structured signal.

Output: results/benchmark/v33/cds_negative_controls.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "cds_tool"))
from rerun_v33 import SEED, load_cohort  # noqa: E402
from cds import circularity_detection_score  # noqa: E402

OUT = REPO / "results" / "benchmark" / "v33" / "cds_negative_controls.json"

FAMILIES = {
    "A: ECM / stroma (no EMT-TF)": {
        "A1 fibrillar collagens": ["COL1A1", "COL1A2", "COL3A1", "COL5A1",
                                   "COL6A1", "COL6A2", "COL6A3", "COL12A1"],
        "A2 fibroblast activation": ["FAP", "PDGFRA", "PDGFRB", "ACTA2",
                                     "THY1", "MME", "CSPG4", "NGFR"],
        "A3 ECM remodeling": ["MMP2", "MMP9", "MMP14", "MMP11", "TIMP1",
                              "TIMP2", "TIMP3", "LOX"],
        "A4 endothelium": ["PECAM1", "VWF", "CDH5", "KDR", "FLT1", "CLDN5",
                           "ESAM", "AQP1"],
        "A5 pericyte": ["RGS5", "NOTCH3", "PDGFRB", "MCAM", "CSPG4",
                        "ACTA2", "TAGLN", "MYL9"],
        "A6 proteoglycan matrix": ["ACAN", "LUM", "DCN", "BGN", "VCAN",
                                   "FMOD", "PRELP", "ASPN"],
        "A7 basement membrane": ["LAMA2", "LAMB1", "LAMC1", "COL4A1",
                                 "COL4A2", "NID1", "NID2", "HSPG2"],
        "A8 myofibroblastic CAF": ["COL11A1", "CTHRC1", "MMP10", "INHBA",
                                   "TGFBI", "TAGLN", "MYL9", "ACTA2"],
        "A9 lymphatic endothelium": ["PROX1", "LYVE1", "PDPN", "CCL21",
                                     "FLT4", "MMRN1", "PDPN", "VEGFC"],
        "A10 protease stroma": ["MMP1", "MMP12", "MMP19", "PLAU", "PLAUR",
                                "SERPINE1", "CTSL", "CTSB"],
    },
    "B: tumor purity / proliferation": {
        "B1 cell cycle core": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM4",
                               "MCM6", "BIRC5", "CCNB1"],
        "B2 E2F targets": ["E2F1", "MCM3", "MCM5", "CDC6", "CDC45", "RRM1",
                           "TYMS", "DHFR"],
        "B3 MYC targets": ["MYC", "MAX", "NPM1", "ODC1", "CDC25A", "SKP2",
                           "EIF4A1", "NCL"],
        "B4 DNA replication/repair": ["PCNA", "FEN1", "LIG1", "POLE2",
                                      "PARP1", "XRCC1", "MRE11", "RAD51"],
        "B5 spindle checkpoint": ["BUB1", "BUB1B", "AURKA", "AURKB", "PLK1",
                                  "CENPF", "TTK", "NDC80"],
        "B6 translation machinery": ["RPS6", "RPL13A", "RPLP0", "EEF1A1",
                                     "EIF3L", "RPS3", "RPL7", "EEF2"],
        "B7 glycolysis": ["SLC2A1", "HK2", "PFKP", "ALDOA", "ENO1", "PKM",
                          "LDHA", "PDK1"],
        "B8 oxidative phosphorylation": ["ATP5F1A", "NDUFS1", "COX4I1",
                                         "UQCRC2", "SDHB", "IDH3A", "CS",
                                         "MDH2"],
        "B9 stress response": ["CDKN1A", "GADD45A", "DDIT3", "ATF3", "JUN",
                               "FOS", "HSPA1A", "DNAJB1"],
        "B10 unfolded protein response": ["HSPA5", "DDIT3", "ATF4", "XBP1",
                                          "PDIA4", "HERPUD1", "ATF6B",
                                          "DNAJB9"],
    },
    "C: lineage / signaling programs": {
        "C1 melanocytic differentiation": ["MLANA", "PMEL", "TYR", "TYRP1",
                                           "DCT", "OCA2", "SLC45A2", "RAB38"],
        "C2 keratinocyte basal": ["KRT5", "KRT14", "KRT17", "KRT16",
                                  "SPRR1B", "SPRR3", "IVL", "KRT1"],
        "C3 hypoxia": ["VEGFA", "CA9", "LDHA", "BNIP3", "PDK1", "SLC2A1",
                       "NDRG1", "DDIT4"],
        "C4 interferon response": ["IFI27", "IFI44L", "IFIT1", "ISG15",
                                   "MX1", "OAS1", "OAS2", "RSAD2"],
        "C5 TNF/NF-kB response": ["NFKBIA", "TNFAIP3", "BIRC3", "TRAF1",
                                  "CD83", "CXCL3", "CXCL2", "CCL2"],
        "C6 antigen presentation I": ["HLA-A", "HLA-B", "HLA-C", "B2M",
                                      "TAP1", "TAP2", "PSMB8", "CD74"],
        "C7 complement": ["C1QA", "C1QB", "C1QC", "C3", "CFB", "C4A",
                          "C1R", "C1S"],
        "C8 senescence": ["CDKN2A", "CDKN1A", "GLB1", "SERPINE1", "IL6",
                          "CXCL8", "MMP3", "MMP1"],
        "C9 Wnt targets": ["WNT5A", "FZD7", "LRP6", "SFRP1", "DKK1",
                           "LEF1", "AXIN2", "TCF7L2"],
        "C10 MHC class II": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DMA",
                             "CD74", "CIITA", "HLA-DOB", "HLA-DMB"],
    },
}

COHORT_JOBS = [("Gide_2019_cBio", "Gide 2019 RECIST (n=73)"),
               ("Hugo_2016", "Hugo 2016 RECIST (n=28)")]


def dense(adata):
    return (adata.X.toarray() if hasattr(adata.X, "toarray")
            else np.asarray(adata.X))


def median_split(score):
    return (score > np.median(score)).astype(int)


def run_cds(X, y, genes, gene_list, rng_note=""):
    r = circularity_detection_score(X, y, genes, gene_list)
    r = json.loads(json.dumps(r, default=float))
    return {"CDS": r["CDS"],
            "C1": r["components"]["C1_MI_endpoint_genes"],
            "C2": r["components"]["C2_rho_endpoint_response"],
            "C3": r["components"]["C3_MI_skew_ratio"],
            "risk": r["risk_level"].split("(")[0].strip(),
            "n_genes_found": len([g for g in genes if g in gene_list]),
            "n_genes_declared": len(genes)}


def main():
    out = {"meta": {
        "seed": SEED,
        "design": "30 gene-program pseudo-endpoints (median split of "
                  "program mean; declared genes = program genes), "
                  "cytolytic positive control (GZMA/PRF1 mean, median "
                  "split), 10 random-label nulls per cohort; cds v1.1.0",
        "note": "pseudo-endpoints are gene-defined by construction, so "
                "HIGH labels indicate program-agnostic detection of "
                "gene-defined surrogates; the random-label null "
                "calibrates the CDS floor",
    }}
    for cname, clabel in COHORT_JOBS:
        adata, _, y_true = load_cohort(cname)
        X = dense(adata)
        gl = list(adata.var_names)
        res = {"pseudo_endpoints": {}, "references": {}}

        for fam, sets in FAMILIES.items():
            for sname, gl_set in sets.items():
                found = [g for g in gl_set if g in gl]
                if len(found) < 3:
                    res["pseudo_endpoints"][f"{fam}|{sname}"] = {
                        "skipped": f"only {len(found)} genes found"}
                    continue
                score = X[:, [gl.index(g) for g in found]].mean(axis=1)
                yp = median_split(score)
                r = run_cds(X, yp, found, gl)
                r["family"] = fam.split(":")[0]
                res["pseudo_endpoints"][f"{fam}|{sname}"] = r
                print(f"  {clabel} {sname}: CDS={r['CDS']:.3f} "
                      f"({r['risk']})", flush=True)

        # positive control: cytolytic surrogate (gene-defined, as in Riaz)
        cyt_genes = [g for g in ["GZMA", "PRF1"] if g in gl]
        cyt = X[:, [gl.index(g) for g in cyt_genes]].mean(axis=1)
        res["references"]["cytolytic_surrogate"] = run_cds(
            X, median_split(cyt), cyt_genes, gl)

        # random-label null
        rng = np.random.RandomState(SEED)
        nulls = []
        for i in range(10):
            yp = rng.permutation(y_true)
            nulls.append(run_cds(X, yp, cyt_genes, gl)["CDS"])
        res["references"]["random_null"] = {
            "n": 10, "cds_values": [round(v, 3) for v in nulls],
            "mean": round(float(np.mean(nulls)), 3),
            "max": round(float(np.max(nulls)), 3)}
        out[clabel] = res

        vals = [v["CDS"] for v in res["pseudo_endpoints"].values()
                if "CDS" in v]
        risks = [v["risk"] for v in res["pseudo_endpoints"].values()
                 if "CDS" in v]
        print(f"[{clabel}] pseudo-endpoints: n={len(vals)} "
              f"CDS mean={np.mean(vals):.3f} min={np.min(vals):.3f} "
              f"max={np.max(vals):.3f} HIGH={risks.count('HIGH')}/"
              f"{len(risks)}", flush=True)

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
