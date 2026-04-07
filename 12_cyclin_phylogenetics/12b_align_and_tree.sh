#!/bin/bash
# 02_align_and_tree.sh
#
# MAFFT alignment, trimAl trimming, IQ-TREE2 phylogeny
# GmJAG1 Manuscript - Task 02
#
# Usage: bash scripts/02_align_and_tree.sh
#   (run from cyclin_phylogeny/ directory)
# Requires: MAFFT, trimAl, IQ-TREE2 loaded via module

set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SEQ_DIR="${BASE_DIR}/01_sequences"
ALN_DIR="${BASE_DIR}/02_alignment"
TREE_DIR="${BASE_DIR}/03_tree"

INPUT="${SEQ_DIR}/CYCD_all_species.fa"
ALIGNED="${ALN_DIR}/CYCD_aligned.fa"
TRIMMED="${ALN_DIR}/CYCD_aligned_trimmed.fa"
TREE_PREFIX="${TREE_DIR}/CYCD_tree"

# Detect iqtree executable name
if command -v iqtree2 &>/dev/null; then
    IQTREE="iqtree2"
elif command -v iqtree &>/dev/null; then
    IQTREE="iqtree"
else
    echo "ERROR: Neither iqtree2 nor iqtree found in PATH"
    exit 1
fi

echo "=========================================="
echo "Alignment & Phylogeny Pipeline"
echo "=========================================="
echo "Input:  ${INPUT}"
echo "IQ-TREE: ${IQTREE} ($(${IQTREE} --version 2>&1 | head -1))"
echo ""

# Verify input exists
if [ ! -f "${INPUT}" ]; then
    echo "ERROR: Input FASTA not found: ${INPUT}"
    echo "Run 01_extract_CYCD.py first."
    exit 1
fi

N_SEQ=$(grep -c '^>' "${INPUT}")
echo "Input sequences: ${N_SEQ}"
echo ""

# ----------------------------------------------------------
# Step 1: MAFFT multiple sequence alignment
# ----------------------------------------------------------
echo "[Step 1] MAFFT alignment..."
mafft --auto --thread 4 --reorder "${INPUT}" > "${ALIGNED}"

ALN_LEN=$(awk '/^>/{next} {print length($0); exit}' "${ALIGNED}")
echo "  Output: ${ALIGNED}"
echo "  Alignment length: ${ALN_LEN} positions"
echo ""

# ----------------------------------------------------------
# Step 2: trimAl - remove poorly aligned regions
# ----------------------------------------------------------
echo "[Step 2] trimAl trimming..."
trimal \
    -in "${ALIGNED}" \
    -out "${TRIMMED}" \
    -automated1 \
    -htmlout "${ALN_DIR}/CYCD_alignment_report.html" \
    -colnumbering > "${ALN_DIR}/CYCD_kept_columns.txt" 2>&1 || true

TRIM_LEN=$(awk '/^>/{next} {print length($0); exit}' "${TRIMMED}")
echo "  Output: ${TRIMMED}"
echo "  Trimmed length: ${TRIM_LEN} positions (${ALN_LEN} -> ${TRIM_LEN})"
echo "  HTML report: ${ALN_DIR}/CYCD_alignment_report.html"
echo ""

# ----------------------------------------------------------
# Step 3: IQ-TREE2 phylogeny
# ----------------------------------------------------------
echo "[Step 3] IQ-TREE2 phylogenetic analysis..."
echo "  Model selection: automatic (TEST)"
echo "  Bootstrap: 1000 ultrafast + 1000 SH-aLRT"
echo ""

${IQTREE} \
    -s "${TRIMMED}" \
    --prefix "${TREE_PREFIX}" \
    -m TEST \
    -bb 1000 \
    -alrt 1000 \
    -nt 4 \
    --seed 42 \
    -redo

echo ""
echo "=========================================="
echo "RESULTS"
echo "=========================================="
echo "  Alignment:       ${ALIGNED}"
echo "  Trimmed:         ${TRIMMED}"
echo "  Tree:            ${TREE_PREFIX}.treefile"
echo "  Best model:      $(grep 'Best-fit model' ${TREE_PREFIX}.iqtree | head -1 || echo 'see .iqtree file')"
echo "  Log-likelihood:  $(grep 'Log-likelihood' ${TREE_PREFIX}.iqtree | head -1 || echo 'see .iqtree file')"
echo ""
echo "Next: run 03_visualize_tree.R"
echo "=========================================="
