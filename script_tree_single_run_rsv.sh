#!/bin/bash
source config.sh

# Function to display step completion
step_complete() {
    echo ""
    echo "========================================"
    echo "Step $1 complete: $2"
    echo "========================================"
    echo ""
}

########################################
# Step 0: Setup and dependency checks
########################################
echo "Step 0: Setup and dependency checks..."

# Make directories
mkdir -p trees 
mkdir -p trees/trees_logs

# Check required tools
command -v mafft >/dev/null 2>&1 || { echo >&2 "Error: MAFFT not found"; exit 1; }
command -v iqtree >/dev/null 2>&1 || { echo >&2 "Error: IQ-TREE not found"; exit 1; }

# Verify input files exist
[ -f "consensus_sequences/RSVA_consensus.fasta" ] || { echo >&2 "Error: RSVA_consensus.fasta not found"; exit 1; }
[ -f "consensus_sequences/RSVB_consensus.fasta" ] || { echo >&2 "Error: RSVB_consensus.fasta not found"; exit 1; }

step_complete "0" "Success!!"

########################################
# Step 5: Small scale MAFFT Alignment with mafft
########################################
echo "Step 5: Small scale MAFFT Alignment with mafft..."

# RSVA Alignment
mafft --auto --thread $THREADS \
      --reorder \
      consensus_sequences/RSVA_consensus.fasta > trees/MAFFT_RSVA.fasta 2> trees/trees_logs/mafft_rsva.log

# RSVB Alignment
mafft --auto --thread $THREADS \
      --reorder \
      consensus_sequences/RSVB_consensus.fasta > trees/MAFFT_RSVB.fasta 2> trees/trees_logs/mafft_rsvb.log

MAFFT_RSVA_COUNT=$(grep -c "^>" trees/MAFFT_RSVA.fasta)
MAFFT_RSVB_COUNT=$(grep -c "^>" trees/MAFFT_RSVB.fasta)

step_complete "5" "Success!! Aligned RSVA: $MAFFT_RSVA_COUNT sequences!! Aligned RSVB: $MAFFT_RSVB_COUNT sequences!!"

########################################
# Step 6: Single-run Tree Construction with iqtree
########################################
echo "Step 6: Single-run Tree Construction with iqtree..."

# RSVA Tree
iqtree -s trees/MAFFT_RSVA.fasta \
       -m MFP -T $THREADS \
       -bb 1000 -alrt 1000 \
       -pre trees/RSVA_tree \
       2> trees/trees_logs/iqtree_mafft_rsva.log

# RSVB Tree
iqtree -s trees/MAFFT_RSVB.fasta \
       -m MFP -T $THREADS \
       -bb 1000 -alrt 1000 \
       -pre trees/RSVB_tree \
       2> trees/trees_logs/iqtree_mafft_rsvb.log

step_complete "6" "Success!!"

echo "========================================"
echo "Single-run tree analysis successfully completed!"
echo "Output files:"
echo "- RSVA Alignment: trees/MAFFT_RSVA.fasta"
echo "- RSVB Alignment: trees/MAFFT_RSVB.fasta"
echo "- RSVA Tree: trees/TREE_MAFFT_RSVA.treefile"
echo "- RSVB Tree: trees/TREE_MAFFT_RSVB.treefile"
echo "========================================"
