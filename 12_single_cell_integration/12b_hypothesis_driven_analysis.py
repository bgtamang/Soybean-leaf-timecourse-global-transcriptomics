#!/usr/bin/env python3
# Script 43: Hypothesis-Driven Single-Cell Analysis
# Tests each link of the JAG1 -> EAR -> TPR -> HDAC repression pathway
# using Fan et al. 2025 soybean shoot apex snRNA-seq data
# Run on HPC: sbatch run_hypothesis.sh (see below for conda env)

import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")
from scipy import sparse
from scipy.stats import mannwhitneyu, spearmanr
from datetime import datetime

# ===== SETUP =====

BASE_DIR = "/home/a-m/bgtamang/Soybean_Leafshape_RNASeq/soybean_RNASeq_2025/data/Fang_et_al"
DATA_DIR = os.path.join(BASE_DIR, "raw_data")
OUT_DIR  = os.path.join(BASE_DIR, "results", "hypothesis_driven")
os.makedirs(OUT_DIR, exist_ok=True)

# HPC environment: conda activate scsoy (Python 3.10, scanpy, scipy)

print("=" * 64)
print("  SCRIPT 43: HYPOTHESIS-DRIVEN SINGLE-CELL ANALYSIS")
print("  JAG1 -> TPR -> HDAC -> Cyclin repression pathway")
print("=" * 64)
print(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 64)
print()


# =====================================================================
# GENE LISTS
# All ZH13 IDs verified against SoyBase pangene lookup table
# (gene_id_conversion_full_lookup.csv). Hyphens match h5ad var_names.
# =====================================================================

# --- Link 1: The Regulator ---
jag_genes = [
    {"name": "GmJAG1", "wm82": "Glyma.20G116200", "zh13": "SoyZH13-20G103500"},
    {"name": "GmJAG2", "wm82": "Glyma.10G273800", "zh13": "SoyZH13-10G253200"},
]

# --- Link 2: The Repression Machinery ---
# 12 TPR co-repressors (manuscript Table S8)
tpr_genes = [
    {"name": "GmTPR1",  "wm82": "Glyma.08G214600", "zh13": "SoyZH13-08G204600", "at_ortholog": "TPL"},
    {"name": "GmTPR2",  "wm82": "Glyma.07G028000", "zh13": "SoyZH13-07G026900", "at_ortholog": "TPL"},
    {"name": "GmTPR3",  "wm82": "Glyma.13G367300", "zh13": "SoyZH13-13G338800", "at_ortholog": "TPR1"},
    {"name": "GmTPR4",  "wm82": "Glyma.15G006100", "zh13": "SoyZH13-15G006300", "at_ortholog": "TPR1"},
    {"name": "GmTPR5",  "wm82": "Glyma.03G233400", "zh13": "SoyZH13-03G212900", "at_ortholog": "TPR2"},
    {"name": "GmTPR6",  "wm82": "Glyma.19G230500", "zh13": "SoyZH13-19G216300", "at_ortholog": "TPR3"},
    {"name": "GmTPR7",  "wm82": "Glyma.20G238400", "zh13": "SoyZH13-20G222000", "at_ortholog": "TPR4"},
    {"name": "GmTPR8",  "wm82": "Glyma.10G150000", "zh13": "SoyZH13-10G136300", "at_ortholog": "TPR3"},
    {"name": "GmTPR9",  "wm82": "Glyma.17G112500", "zh13": "SoyZH13-17G108500", "at_ortholog": "TPL family"},
    {"name": "GmTPR10", "wm82": "Glyma.13G158800", "zh13": "SoyZH13-13G141600", "at_ortholog": "TPR1"},
    {"name": "GmTPR11", "wm82": "Glyma.04G065000", "zh13": "SoyZH13-04G061700", "at_ortholog": "Soybean-specific"},
    {"name": "GmTPR12", "wm82": "Glyma.06G066200", "zh13": "SoyZH13-06G061800", "at_ortholog": "Soybean-specific"},
]

# 15 HDAC complex genes (curated from Yang et al. 2018; 9 catalytic + 6 subunits).
# Removed 5 misannotated genes from original 20-gene keyword-grep list:
#   Glyma.01G243000 (IST1P/UBP), Glyma.05G185800 (Sm D3/SNRPD3),
#   Glyma.08G143900 (Sm D3/SNRPD3), Glyma.12G203000 (IST1P/ESCRT-III),
#   Glyma.13G298700 (IST1P/ESCRT-III)
hda_genes = [
    {"name": "GmHDA1",  "wm82": "Glyma.03G019600", "zh13": None,                    "at_ortholog": "SAP18"},
    {"name": "GmHDA2",  "wm82": "Glyma.04G033900", "zh13": None,                    "at_ortholog": "SAP30L"},
    {"name": "GmHDA3",  "wm82": "Glyma.04G034000", "zh13": "SoyZH13-04G032000",     "at_ortholog": "SAP30L"},
    {"name": "GmHDA4",  "wm82": "Glyma.05G012850", "zh13": None,                    "at_ortholog": "HDA15"},
    {"name": "GmHDA5",  "wm82": "Glyma.05G021400", "zh13": "SoyZH13-05G020400",     "at_ortholog": "RPD3"},
    {"name": "GmHDA6",  "wm82": "Glyma.05G192600", "zh13": "SoyZH13-05G180300",     "at_ortholog": "HDA6"},
    {"name": "GmHDA7",  "wm82": "Glyma.06G034100", "zh13": "SoyZH13-06G032300",     "at_ortholog": "SAP30L"},
    {"name": "GmHDA8",  "wm82": "Glyma.06G034200", "zh13": "SoyZH13-06G032400",     "at_ortholog": "SAP30L"},
    {"name": "GmHDA9",  "wm82": "Glyma.07G081500", "zh13": "SoyZH13-07G076000",     "at_ortholog": "SAP18"},
    {"name": "GmHDA10", "wm82": "Glyma.11G187800", "zh13": "SoyZH13-11G159700",     "at_ortholog": "HDAC1"},
    {"name": "GmHDA11", "wm82": "Glyma.12G086700", "zh13": "SoyZH13-12G079300",     "at_ortholog": "HDAC1"},
    {"name": "GmHDA12", "wm82": "Glyma.12G188200", "zh13": "SoyZH13-12G168600",     "at_ortholog": "Predicted HDA"},
    {"name": "GmHDA13", "wm82": "Glyma.17G078000", "zh13": "SoyZH13-17G075400",     "at_ortholog": "RPD3"},
    {"name": "GmHDA14", "wm82": "Glyma.17G120900", "zh13": "SoyZH13-17G116900",     "at_ortholog": "HDA15"},
    {"name": "GmHDA15", "wm82": "Glyma.17G229600", "zh13": "SoyZH13-17G218400",     "at_ortholog": "HDAC11"},
]

# --- Link 3: KRPs — NOT JAG1 targets in soybean ---
krp_genes = [
    {"name": "GmKRP1a", "wm82": "Glyma.01G185400", "zh13": "SoyZH13-01G168700"},
    {"name": "GmKRP1b", "wm82": "Glyma.11G056700", "zh13": "SoyZH13-11G055900"},
    {"name": "GmKRP2a", "wm82": "Glyma.05G085000", "zh13": "SoyZH13-05G084200"},
    {"name": "GmKRP2b", "wm82": "Glyma.17G175700", "zh13": "SoyZH13-17G167700"},
    {"name": "GmKRP3",  "wm82": "Glyma.08G354300", "zh13": "SoyZH13-08G333700"},
    {"name": "GmKRP4",  "wm82": "Glyma.20G198800", "zh13": "SoyZH13-20G183900"},
    {"name": "GmKRP5",  "wm82": "Glyma.18G170800", "zh13": "SoyZH13-18G147800"},
    {"name": "GmKRP6",  "wm82": "Glyma.02G133700", "zh13": "SoyZH13-02G124700"},
    {"name": "GmKRP7",  "wm82": "Glyma.07G211300", "zh13": "SoyZH13-07G196800"},
]

# --- Link 4: JAG1-target cyclins (7, all de-repressed in narrow leaf) ---
cyclin_targets = [
    {"name": "CYCLIN-U4-1a", "wm82": "Glyma.03G107800", "zh13": "SoyZH13-03G096000", "type": "U",  "bulk_logFC": 0.94},
    {"name": "CYCLIN-A3-1",  "wm82": "Glyma.06G044000", "zh13": "SoyZH13-06G041000", "type": "A",  "bulk_logFC": 1.74},
    {"name": "CYCLIN-U4-1b", "wm82": "Glyma.07G119100", "zh13": "SoyZH13-07G112200", "type": "U",  "bulk_logFC": 1.50},
    {"name": "CYCLIN-D1-1",  "wm82": "Glyma.08G291000", "zh13": "SoyZH13-08G274600", "type": "D",  "bulk_logFC": 1.16},
    {"name": "CYCD6",        "wm82": "Glyma.13G247600", "zh13": "SoyZH13-13G224900", "type": "D",  "bulk_logFC": 1.65},
    {"name": "CYCLIN-U1-1",  "wm82": "Glyma.13G306600", "zh13": "SoyZH13-13G280200", "type": "U",  "bulk_logFC": 0.60},
    {"name": "CYCLIN-H",     "wm82": "Glyma.17G201300", "zh13": "SoyZH13-17G189900", "type": "H",  "bulk_logFC": 0.68},
]

# --- Link 5: Core CDKs — constitutive machinery, not JAG1-regulated ---
cdk_genes = [
    {"name": "CDKA-1",  "wm82": "Glyma.05G123500", "zh13": None,                    "class": "CDKA"},
    {"name": "CDKA-2",  "wm82": "Glyma.09G030000", "zh13": "SoyZH13-09G029200", "class": "CDKA"},
    {"name": "CDKA-3",  "wm82": "Glyma.08G078600", "zh13": "SoyZH13-08G074700", "class": "CDKA"},
    {"name": "CDKA-4",  "wm82": "Glyma.15G135800", "zh13": "SoyZH13-15G131000", "class": "CDKA"},
    {"name": "CDKB1",   "wm82": "Glyma.07G021100", "zh13": "SoyZH13-07G019900", "class": "CDKB1"},
    {"name": "CDKB2-1", "wm82": "Glyma.14G224000", "zh13": "SoyZH13-14G206300", "class": "CDKB2"},
    {"name": "CDKB2-2", "wm82": "Glyma.17G262300", "zh13": "SoyZH13-17G249500", "class": "CDKB2"},
    {"name": "CDKB2-3", "wm82": "Glyma.09G071300", "zh13": "SoyZH13-09G067900", "class": "CDKB2"},
    {"name": "CDKB2-4", "wm82": "Glyma.07G070200", "zh13": "SoyZH13-07G065700", "class": "CDKB2"},
]


# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

def get_gene_expression(adata, gene_id):
    """Extract expression vector for a single gene from the AnnData object."""
    if gene_id is None or gene_id not in adata.var_names:
        return None
    idx = list(adata.var_names).index(gene_id)
    if sparse.issparse(adata.X):
        return np.array(adata.X[:, idx].toarray()).flatten()
    return np.array(adata.X[:, idx]).flatten()


def cluster_stats(adata, gene_id, cluster_col):
    """Per-cluster expression stats for one gene."""
    expr = get_gene_expression(adata, gene_id)
    if expr is None:
        return None
    rows = []
    for ct in sorted(adata.obs[cluster_col].unique()):
        mask = adata.obs[cluster_col].values == ct
        ct_expr = expr[mask]
        n = len(ct_expr)
        n_pos = int(np.sum(ct_expr > 0))
        rows.append({
            "cluster": ct, "n_cells": n, "n_expressing": n_pos,
            "pct_expressing": round(100 * n_pos / n, 2) if n > 0 else 0,
            "mean_all_cells": round(float(np.mean(ct_expr)), 4),
        })
    return pd.DataFrame(rows)


def gene_family_per_cluster(adata, gene_list, cluster_col, family_name):
    """Map a gene family across all SAM clusters. Prints per-gene summary."""
    all_results = []
    found, skipped, not_found = 0, 0, 0

    for gene in gene_list:
        gid = gene.get("zh13")
        name = gene["name"]
        if gid is None:
            print(f"    {name} ({gene['wm82']}): No ZH13 ID -- SKIPPED")
            skipped += 1
            continue
        df = cluster_stats(adata, gid, cluster_col)
        if df is None:
            print(f"    {name} ({gid}): NOT FOUND in matrix")
            not_found += 1
            continue

        df["gene_name"] = name
        df["gene_id_zh13"] = gid
        df["gene_id_wm82"] = gene["wm82"]
        for extra in ["at_ortholog", "type", "class", "bulk_logFC"]:
            if extra in gene:
                df[extra] = gene[extra]
        all_results.append(df)
        found += 1

        # One-line summary
        total_pct = 100 * df["n_expressing"].sum() / df["n_cells"].sum()
        top = df.loc[df["pct_expressing"].idxmax()]
        lp1 = df[df["cluster"] == "Leaf primordium1"]
        lp1_pct = lp1["pct_expressing"].values[0] if len(lp1) > 0 else 0
        tag = ""
        if "at_ortholog" in gene:
            tag = f", {gene['at_ortholog']}"
        elif "type" in gene:
            tag = f", {gene['type']}-type"
        elif "class" in gene:
            tag = f", {gene['class']}"
        print(f"    {name}{tag}: {total_pct:.1f}% overall, LP1={lp1_pct:.1f}%, "
              f"top={top['cluster']} ({top['pct_expressing']:.1f}%)")

    print(f"  [{family_name}: {found} found, {skipped} skipped (no ZH13), {not_found} not in matrix]")
    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return None


def mean_pct_per_cluster(adata, gene_list, cluster_col):
    """Mean % expressing across a gene family for each cluster."""
    clusters = sorted(adata.obs[cluster_col].unique())
    pcts = {cl: [] for cl in clusters}
    for gene in gene_list:
        gid = gene.get("zh13")
        if gid is None:
            continue
        expr = get_gene_expression(adata, gid)
        if expr is None:
            continue
        for cl in clusters:
            mask = adata.obs[cluster_col].values == cl
            n = mask.sum()
            if n > 0:
                pcts[cl].append(100 * np.sum(expr[mask] > 0) / n)
    return {cl: np.mean(vals) if vals else 0 for cl, vals in pcts.items()}


# ===== SECTION 1: LOAD DATA =====

print("========================================")
print("SECTION 1: LOAD DATA")
print("========================================")
print()

sam = sc.read_h5ad(os.path.join(DATA_DIR, "SAM.scRNA.h5ad"))
leaf = sc.read_h5ad(os.path.join(DATA_DIR, "Leaf.scRNA.h5ad"))
cluster_col = "celltype"

print(f"  SAM:  {sam.n_obs} cells, {sam.n_vars} genes, {len(sam.obs[cluster_col].unique())} clusters")
print(f"  Leaf: {leaf.n_obs} cells, {leaf.n_vars} genes, {len(leaf.obs[cluster_col].unique())} clusters")
clusters = sorted(sam.obs[cluster_col].unique())
print(f"  SAM clusters: {clusters}")
print()


# ===== SECTION 2: LINK 1 — JAG1/JAG2 =====

print("========================================")
print("SECTION 2: LINK 1 -- JAG1 / JAG2")
print("========================================")
print()

jag_df = gene_family_per_cluster(sam, jag_genes, cluster_col, "JAG")
if jag_df is not None:
    jag_df.to_csv(os.path.join(OUT_DIR, "link1_jag_per_cluster.csv"), index=False)

# SAM vs Leaf comparison
print("\n  JAG1/JAG2 tissue comparison:")
for gene in jag_genes:
    for adata, label in [(sam, "SAM"), (leaf, "Leaf")]:
        expr = get_gene_expression(adata, gene["zh13"])
        if expr is not None:
            pct = 100 * np.sum(expr > 0) / len(expr)
            print(f"    {gene['name']} in {label}: {pct:.2f}% ({int(np.sum(expr > 0))}/{len(expr)} cells)")
print()


# ===== SECTION 3: LINK 2 — TPR + HDAC =====

print("========================================")
print("SECTION 3: LINK 2 -- TPR + HDAC")
print("========================================")
print()

print("  --- TPR co-repressors (12 genes) ---")
tpr_df = gene_family_per_cluster(sam, tpr_genes, cluster_col, "TPR")
if tpr_df is not None:
    tpr_df.to_csv(os.path.join(OUT_DIR, "link2_tpr_per_cluster.csv"), index=False)

print("\n  --- HDAC complex (15 genes, curated) ---")
hda_df = gene_family_per_cluster(sam, hda_genes, cluster_col, "HDAC")
if hda_df is not None:
    hda_df.to_csv(os.path.join(OUT_DIR, "link2_hda_per_cluster.csv"), index=False)
print()


# ===== SECTION 4: LINK 3 — KRPs (NEGATIVE CONTROL) =====

print("========================================")
print("SECTION 4: LINK 3 -- KRPs (negative control)")
print("========================================")
print()

print("  --- KRP expression in SAM (9 genes) ---")
krp_df = gene_family_per_cluster(sam, krp_genes, cluster_col, "KRP")
if krp_df is not None:
    krp_df.to_csv(os.path.join(OUT_DIR, "link3_krp_per_cluster.csv"), index=False)

# KRP SAM vs Leaf — tests endoreduplication signal
print("\n  --- KRP SAM vs Leaf comparison ---")
krp_tissue = []
for gene in krp_genes:
    gid = gene["zh13"]
    sam_expr = get_gene_expression(sam, gid)
    leaf_expr = get_gene_expression(leaf, gid)
    if sam_expr is not None and leaf_expr is not None:
        sam_pct = 100 * np.sum(sam_expr > 0) / len(sam_expr)
        leaf_pct = 100 * np.sum(leaf_expr > 0) / len(leaf_expr)
        ratio = leaf_pct / sam_pct if sam_pct > 0 else float('inf')
        direction = "HIGHER in Leaf" if ratio > 1.5 else ("LOWER in Leaf" if ratio < 0.67 else "SIMILAR")
        stat, pval = mannwhitneyu(leaf_expr, sam_expr, alternative='two-sided')
        print(f"    {gene['name']}: SAM={sam_pct:.2f}%, Leaf={leaf_pct:.2f}%, "
              f"ratio={ratio:.2f}, p={pval:.2e} -> {direction}")
        krp_tissue.append({
            "gene": gene["name"], "sam_pct": sam_pct, "leaf_pct": leaf_pct,
            "ratio": round(ratio, 2), "pvalue": pval, "direction": direction
        })

krp_tissue_df = pd.DataFrame(krp_tissue)
krp_tissue_df.to_csv(os.path.join(OUT_DIR, "link3_krp_sam_vs_leaf.csv"), index=False)
n_higher = sum(1 for r in krp_tissue if r["direction"] == "HIGHER in Leaf")
print(f"\n  Result: {n_higher}/{len(krp_tissue)} KRPs higher in Leaf")
print()


# ===== SECTION 5: LINK 4 — CYCLINS (KEY TEST) =====

print("========================================")
print("SECTION 5: LINK 4 -- JAG1-target cyclins")
print("========================================")
print("  7 cyclins, all upregulated in narrow (de-repressed when JAG1 non-functional)")
print()

cyc_df = gene_family_per_cluster(sam, cyclin_targets, cluster_col, "Cyclin targets")
if cyc_df is not None:
    cyc_df.to_csv(os.path.join(OUT_DIR, "link4_cyclin_targets_per_cluster.csv"), index=False)
print()


# ===== SECTION 6: LINK 5 — CDKs (SECOND NEGATIVE CONTROL) =====

print("========================================")
print("SECTION 6: LINK 5 -- CDKs (constitutive)")
print("========================================")
print()

cdk_df = gene_family_per_cluster(sam, cdk_genes, cluster_col, "CDK")
if cdk_df is not None:
    cdk_df.to_csv(os.path.join(OUT_DIR, "link5_cdk_per_cluster.csv"), index=False)
print()


# ===== SECTION 7: SUMMARY TABLE =====

print("========================================")
print("SECTION 7: SUMMARY -- ALL COMPONENTS")
print("========================================")
print()

jag1_expr = get_gene_expression(sam, jag_genes[0]["zh13"])
tpr_pcts = mean_pct_per_cluster(sam, tpr_genes, cluster_col)
hda_pcts = mean_pct_per_cluster(sam, hda_genes, cluster_col)
krp_pcts = mean_pct_per_cluster(sam, krp_genes, cluster_col)
cyc_pcts = mean_pct_per_cluster(sam, cyclin_targets, cluster_col)
cdk_pcts = mean_pct_per_cluster(sam, cdk_genes, cluster_col)

summary_rows = []
for cl in clusters:
    mask = sam.obs[cluster_col].values == cl
    n = int(mask.sum())
    jag1_pct = 100 * np.sum(jag1_expr[mask] > 0) / n if jag1_expr is not None and n > 0 else 0
    summary_rows.append({
        "cluster": cl, "n_cells": n,
        "jag1_pct": round(jag1_pct, 2),
        "tpr_mean_pct": round(tpr_pcts.get(cl, 0), 2),
        "hda_mean_pct": round(hda_pcts.get(cl, 0), 2),
        "krp_mean_pct": round(krp_pcts.get(cl, 0), 2),
        "cyclin_target_mean_pct": round(cyc_pcts.get(cl, 0), 2),
        "cdk_mean_pct": round(cdk_pcts.get(cl, 0), 2),
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUT_DIR, "summary_all_components_by_cluster.csv"), index=False)
print(summary_df.to_string(index=False))


# ===== SECTION 8: INTERPRETATION =====

print()
print("========================================")
print("SECTION 8: INTERPRETATION")
print("========================================")
print()

lp1 = summary_df[summary_df["cluster"] == "Leaf primordium1"].iloc[0]
others = summary_df[summary_df["cluster"] != "Leaf primordium1"]

print(f"  LP1 (n={lp1['n_cells']} cells):")
print(f"    JAG1:     {lp1['jag1_pct']:.1f}%")
print(f"    TPR:      {lp1['tpr_mean_pct']:.1f}%")
print(f"    HDAC:     {lp1['hda_mean_pct']:.1f}%   (machinery)")
print(f"    KRP:      {lp1['krp_mean_pct']:.1f}%")
print(f"    Cyclins:  {lp1['cyclin_target_mean_pct']:.2f}%")
print(f"    CDK:      {lp1['cdk_mean_pct']:.2f}%")

print(f"\n  Other clusters (mean +/- SD, n={len(others)}):")
for col, label in [("jag1_pct", "JAG1"), ("tpr_mean_pct", "TPR"),
                    ("hda_mean_pct", "HDAC"), ("krp_mean_pct", "KRP"),
                    ("cyclin_target_mean_pct", "Cyclins"), ("cdk_mean_pct", "CDK")]:
    m, s = others[col].mean(), others[col].std()
    ratio = lp1[col] / m if m > 0 else float('inf')
    print(f"    {label:10s}: {m:.2f}% +/- {s:.2f}%  (LP1/others = {ratio:.1f}x)")

# LP1 rank per component
print("\n  LP1 rank among 16 clusters:")
for col, label in [("jag1_pct", "JAG1"), ("tpr_mean_pct", "TPR"),
                    ("hda_mean_pct", "HDAC"), ("krp_mean_pct", "KRP"),
                    ("cyclin_target_mean_pct", "Cyclins"), ("cdk_mean_pct", "CDK")]:
    rank = int((summary_df[col].rank(ascending=False))[summary_df["cluster"] == "Leaf primordium1"].values[0])
    print(f"    {label:10s}: #{rank}/16")

# Spearman across 16 clusters
print("\n  Spearman (JAG1 vs each family across 16 clusters):")
jag1_vals = summary_df["jag1_pct"].values
for col, label in [("tpr_mean_pct", "TPR"), ("hda_mean_pct", "HDAC"),
                    ("krp_mean_pct", "KRP"), ("cyclin_target_mean_pct", "Cyclins"),
                    ("cdk_mean_pct", "CDK")]:
    rho, pval = spearmanr(jag1_vals, summary_df[col].values)
    print(f"    JAG1 vs {label:10s}: rho={rho:+.3f}, p={pval:.4f}")


# ===== DONE =====

print()
print("=" * 64)
print("  HYPOTHESIS-DRIVEN ANALYSIS COMPLETE")
print(f"  Results: {OUT_DIR}")
print(f"  Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 64)
