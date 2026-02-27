#!/usr/bin/env python3
# Script 44: Data-Driven Single-Cell Analysis
# Maps all 1,567 JAG1 target genes onto SAM cell types
# Tests whether targets co-localize with JAG1 across clusters
# Requires: 12a output (jag1_targets_1567_converted.csv)

import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from datetime import datetime

# ===== SETUP =====

BASE_DIR = "/home/a-m/bgtamang/Soybean_Leafshape_RNASeq/soybean_RNASeq_2025/data/Fang_et_al"
DATA_DIR = os.path.join(BASE_DIR, "raw_data")
GENE_LIST_DIR = os.path.join(BASE_DIR, "gene_lists")
OUT_DIR = os.path.join(BASE_DIR, "results", "data_driven")
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 64)
print("  SCRIPT 44: DATA-DRIVEN SINGLE-CELL ANALYSIS")
print("  Mapping 1,567 JAG1 targets onto SAM cell types")
print("=" * 64)
print(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 64)
print()


# ===== SECTION 1: LOAD TARGET GENE IDS =====

print("========================================")
print("SECTION 1: LOAD TARGET GENE IDS")
print("========================================")
print()

# Pre-converted by 12a_gene_id_conversion.R
target_file = os.path.join(GENE_LIST_DIR, "jag1_targets_1567_converted.csv")
targets = pd.read_csv(target_file)
print(f"  Loaded {len(targets)} JAG1 target genes")
print(f"  Tiers: {targets['Confidence_Tier'].value_counts().to_dict()}")

# Clean ZH13 IDs: strip whitespace, take first if comma-separated, replace _ with -
targets["GeneID_clean"] = targets["Wm82"].str.strip().str.strip('"')

def clean_zh13(x):
    if pd.isna(x) or str(x).strip() in ('', 'NA', 'nan'):
        return None
    zh13 = str(x).strip().strip('"')
    if "," in zh13:
        zh13 = zh13.split(",")[0].strip()  # multi-mapping: take first
    return zh13.replace("_", "-")  # h5ad uses hyphens

targets["zh13_id"] = targets["Zh13"].apply(clean_zh13)

n_converted = targets["zh13_id"].notna().sum()
n_missing = targets["zh13_id"].isna().sum()
print(f"\n  With ZH13 ID: {n_converted}/{len(targets)} ({100*n_converted/len(targets):.1f}%)")
print(f"  No ZH13 ID:   {n_missing}")

for tier in ["Gold", "Silver", "Bronze"]:
    t = targets[targets["Confidence_Tier"] == tier]
    conv = t["zh13_id"].notna().sum()
    print(f"    {tier}: {conv}/{len(t)} converted")

targets_conv = targets[targets["zh13_id"].notna()].copy()
print(f"\n  Proceeding with {len(targets_conv)} targets")
print()


# ===== SECTION 2: LOAD SAM DATA =====

print("========================================")
print("SECTION 2: LOAD SAM DATA")
print("========================================")
print()

sam = sc.read_h5ad(os.path.join(DATA_DIR, "SAM.scRNA.h5ad"))
cluster_col = "celltype"
clusters = sorted(sam.obs[cluster_col].unique())
print(f"  SAM: {sam.n_obs} cells, {sam.n_vars} genes")
print(f"  Clusters: {clusters}")

# How many targets are in the expression matrix?
in_matrix = targets_conv["zh13_id"].isin(sam.var_names).sum()
print(f"\n  Targets in SAM matrix: {in_matrix}/{len(targets_conv)}")
print(f"  Not in matrix: {len(targets_conv) - in_matrix}")

targets_final = targets_conv[targets_conv["zh13_id"].isin(sam.var_names)].copy()
print(f"  Final target set: {len(targets_final)} genes")
print()


# ===== SECTION 3: TARGET SIGNATURE SCORE PER CLUSTER =====

print("========================================")
print("SECTION 3: TARGET SIGNATURE PER CLUSTER")
print("========================================")
print()

# Get JAG1 expression for reference
jag1_zh13 = "SoyZH13-20G103500"
jag1_expr = None
if jag1_zh13 in sam.var_names:
    idx = list(sam.var_names).index(jag1_zh13)
    if sparse.issparse(sam.X):
        jag1_expr = np.array(sam.X[:, idx].toarray()).flatten()
    else:
        jag1_expr = np.array(sam.X[:, idx]).flatten()

# Extract full target matrix (cells x targets) in one go
target_zh13_ids = targets_final["zh13_id"].tolist()
target_indices = [list(sam.var_names).index(gid) for gid in target_zh13_ids]
print(f"  Extracting expression for {len(target_indices)} targets...")

if sparse.issparse(sam.X):
    target_matrix = np.array(sam.X[:, target_indices].toarray())
else:
    target_matrix = np.array(sam.X[:, target_indices])

# Per-cluster statistics
cluster_scores = []
for cl in clusters:
    mask = sam.obs[cluster_col].values == cl
    n = int(mask.sum())
    cl_matrix = target_matrix[mask, :]

    pct_per_gene = 100 * np.sum(cl_matrix > 0, axis=0) / n
    mean_pct_expressing = float(np.mean(pct_per_gene))
    mean_expr = float(np.mean(cl_matrix))
    n_detected = int(np.sum(np.any(cl_matrix > 0, axis=0)))

    jag1_pct = 0
    if jag1_expr is not None:
        jag1_pct = 100 * np.sum(jag1_expr[mask] > 0) / n if n > 0 else 0

    cluster_scores.append({
        "cluster": cl, "n_cells": n,
        "jag1_pct": round(jag1_pct, 2),
        "n_targets_detected": n_detected,
        "pct_targets_detected": round(100 * n_detected / len(target_indices), 1),
        "target_mean_pct_expressing": round(mean_pct_expressing, 3),
        "target_mean_expression": round(mean_expr, 6),
    })

scores_df = pd.DataFrame(cluster_scores)
scores_df.to_csv(os.path.join(OUT_DIR, "target_score_by_cluster.csv"), index=False)
print()
print(scores_df.to_string(index=False))
print()


# ===== SECTION 4: TIER-SPECIFIC ANALYSIS =====

print("========================================")
print("SECTION 4: TIER-SPECIFIC SCORES")
print("========================================")
print()

tier_results = []
for tier in ["Gold", "Silver", "Bronze"]:
    tier_genes = targets_final[targets_final["Confidence_Tier"] == tier]
    tier_zh13 = tier_genes["zh13_id"].tolist()
    tier_idx = [list(sam.var_names).index(gid) for gid in tier_zh13 if gid in sam.var_names]

    if not tier_idx:
        print(f"  {tier}: No genes in matrix")
        continue

    if sparse.issparse(sam.X):
        tier_matrix = np.array(sam.X[:, tier_idx].toarray())
    else:
        tier_matrix = np.array(sam.X[:, tier_idx])

    print(f"  --- {tier} tier ({len(tier_idx)} genes) ---")
    for cl in clusters:
        mask = sam.obs[cluster_col].values == cl
        n = int(mask.sum())
        cl_mat = tier_matrix[mask, :]
        pct_per_gene = 100 * np.sum(cl_mat > 0, axis=0) / n
        mean_pct = float(np.mean(pct_per_gene))
        mean_expr = float(np.mean(cl_mat))
        jag1_pct = 100 * np.sum(jag1_expr[mask] > 0) / n if jag1_expr is not None and n > 0 else 0

        tier_results.append({
            "tier": tier, "cluster": cl, "n_cells": n, "n_genes": len(tier_idx),
            "jag1_pct": round(jag1_pct, 2),
            "target_mean_pct": round(mean_pct, 3),
            "target_mean_expr": round(mean_expr, 6),
        })

    # LP1 vs others for this tier
    tier_rows = [r for r in tier_results if r["tier"] == tier]
    lp1_r = [r for r in tier_rows if r["cluster"] == "Leaf primordium1"]
    oth_r = [r for r in tier_rows if r["cluster"] != "Leaf primordium1"]
    if lp1_r and oth_r:
        lp1_pct = lp1_r[0]["target_mean_pct"]
        oth_mean = np.mean([r["target_mean_pct"] for r in oth_r])
        if oth_mean > 0:
            print(f"    LP1: {lp1_pct:.3f}% | Others: {oth_mean:.3f}% | Ratio: {lp1_pct/oth_mean:.2f}x")

tier_df = pd.DataFrame(tier_results)
tier_df.to_csv(os.path.join(OUT_DIR, "target_score_by_tier_cluster.csv"), index=False)
print()


# ===== SECTION 5: CORRELATION TESTS =====

print("========================================")
print("SECTION 5: JAG1 vs TARGET CORRELATION")
print("========================================")
print()

jag1_vals = scores_df["jag1_pct"].values
target_pct_vals = scores_df["target_mean_pct_expressing"].values
target_expr_vals = scores_df["target_mean_expression"].values

rho_pct, p_pct = spearmanr(jag1_vals, target_pct_vals)
rho_expr, p_expr = spearmanr(jag1_vals, target_expr_vals)
r_pct, pr_pct = pearsonr(jag1_vals, target_pct_vals)
r_expr, pr_expr = pearsonr(jag1_vals, target_expr_vals)

print(f"  Spearman (JAG1% vs target %):   rho = {rho_pct:+.3f}, p = {p_pct:.4f}")
print(f"  Spearman (JAG1% vs target expr): rho = {rho_expr:+.3f}, p = {p_expr:.4f}")
print(f"  Pearson  (JAG1% vs target %):   r   = {r_pct:+.3f}, p = {pr_pct:.4f}")
print(f"  Pearson  (JAG1% vs target expr): r   = {r_expr:+.3f}, p = {pr_expr:.4f}")

# Per-tier Spearman
print("\n  Tier-specific (Spearman):")
for tier in ["Gold", "Silver", "Bronze"]:
    t_rows = tier_df[tier_df["tier"] == tier]
    if len(t_rows) == len(clusters):
        tier_pct = t_rows.sort_values("cluster")["target_mean_pct"].values
        jag1_sorted = t_rows.sort_values("cluster")["jag1_pct"].values
        rho_t, p_t = spearmanr(jag1_sorted, tier_pct)
        print(f"    {tier:8s}: rho = {rho_t:+.3f}, p = {p_t:.4f}")
print()


# ===== SECTION 6: TOP TARGETS IN LP1 =====

print("========================================")
print("SECTION 6: TOP TARGETS IN LP1")
print("========================================")
print()

lp1_mask = sam.obs[cluster_col].values == "Leaf primordium1"
n_lp1 = int(lp1_mask.sum())
lp1_matrix = target_matrix[lp1_mask, :]

target_lp1_stats = []
for i, (_, row) in enumerate(targets_final.iterrows()):
    col_expr = lp1_matrix[:, i]
    n_pos = int(np.sum(col_expr > 0))
    pct = 100 * n_pos / n_lp1 if n_lp1 > 0 else 0
    target_lp1_stats.append({
        "gene_wm82": row["GeneID_clean"],
        "gene_zh13": row["zh13_id"],
        "tier": row["Confidence_Tier"],
        "pattern": row.get("Pattern", ""),
        "lp1_pct_expressing": round(pct, 2),
        "lp1_mean_expr": round(float(np.mean(col_expr)), 4),
    })

lp1_stats_df = pd.DataFrame(target_lp1_stats)
lp1_stats_df = lp1_stats_df.sort_values("lp1_pct_expressing", ascending=False)
lp1_stats_df.to_csv(os.path.join(OUT_DIR, "all_targets_lp1_expression.csv"), index=False)

print(f"  LP1: {n_lp1} cells")
print(f"  Top 20 targets:")
print(lp1_stats_df.head(20).to_string(index=False))
print()


# ===== SECTION 7: SUMMARY =====

print("========================================")
print("SECTION 7: SUMMARY")
print("========================================")
print()

lp1_score = scores_df[scores_df["cluster"] == "Leaf primordium1"].iloc[0]
others = scores_df[scores_df["cluster"] != "Leaf primordium1"]

print(f"  LP1 target signature:")
print(f"    Targets detected: {lp1_score['n_targets_detected']}/{len(targets_final)} ({lp1_score['pct_targets_detected']}%)")
print(f"    Mean % expressing: {lp1_score['target_mean_pct_expressing']:.3f}%")
print(f"    Mean expression:   {lp1_score['target_mean_expression']:.6f}")

print(f"\n  Other clusters (mean +/- SD):")
print(f"    Mean % expressing: {others['target_mean_pct_expressing'].mean():.3f}% +/- {others['target_mean_pct_expressing'].std():.3f}%")
print(f"    Mean expression:   {others['target_mean_expression'].mean():.6f} +/- {others['target_mean_expression'].std():.6f}")

lp1_ratio = lp1_score['target_mean_pct_expressing'] / others['target_mean_pct_expressing'].mean() if others['target_mean_pct_expressing'].mean() > 0 else 0
print(f"\n  LP1 / Others ratio: {lp1_ratio:.2f}x")

# Targets with zero expression in LP1
n_zero_lp1 = int(np.sum(np.all(lp1_matrix == 0, axis=0)))
print(f"  Targets with zero expression in LP1: {n_zero_lp1}/{len(targets_final)} ({100*n_zero_lp1/len(targets_final):.1f}%)")


# ===== DONE =====

print()
print("=" * 64)
print("  DATA-DRIVEN ANALYSIS COMPLETE")
print(f"  Results: {OUT_DIR}")
print(f"  Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 64)
