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

# Make directories & Add config.sh here
ALL_RUNS="all_runs_trees"
mkdir -p $ALL_RUNS
mkdir -p $ALL_RUNS/trees_logs

# Check for required tools
command -v mafft >/dev/null 2>&1 || { echo >&2 "Error: mafft not found"; exit 1; }
command -v iqtree >/dev/null 2>&1 || { echo >&2 "Error: iqtree not found"; exit 1; }

step_complete "0" "Success!!"

########################################
# Step 5: Combine All Consensus Files Across Runs
########################################
echo "Step 5: Combine All Consensus Files Across Runs..."

# Verify at least some input files exist
if ! compgen -G "*/consensus_sequences/RSVA_consensus.fasta" > /dev/null; then
    echo "ERROR: No RSVA consensus files found" >&2
    exit 1
fi

# RSVA
find . -name "RSVA_consensus.fasta" -exec cat {} + > $ALL_RUNS/ALL_RSVA.fasta
    
# RSVB 
find . -name "RSVB_consensus.fasta" -exec cat {} + > $ALL_RUNS/ALL_RSVB.fasta

# Verify output
RSVA_TOTAL=$(grep -c "^>" $ALL_RUNS/ALL_RSVA.fasta)
RSVB_TOTAL=$(grep -c "^>" $ALL_RUNS/ALL_RSVB.fasta)

step_complete "5" "Success!! Total RSVA: $RSVA_TOTAL sequences!! Total RSVB: $RSVB_TOTAL sequences!!"

########################################
# Step 6: Large-scale MAFFT Alignment with mafft
########################################
echo "Step 6: Large-scale MAFFT Alignment with mafft..."

# RSVA Alignment
mafft --auto --thread $THREADS \
      --large --reorder \
      $ALL_RUNS/ALL_RSVA.fasta > $ALL_RUNS/MAFFT_ALL_RSVA.fasta 2> $ALL_RUNS/trees_logs/mafft_all_rsva.log

# RSVB Alignment
mafft --auto --thread $THREADS \
      --large --reorder \
      $ALL_RUNS/ALL_RSVB.fasta > $ALL_RUNS/MAFFT_ALL_RSVB.fasta 2> $ALL_RUNS/trees_logs/mafft_all_rsvb.log

MAFFT_ALL_RSVA_COUNT=$(grep -c "^>" $ALL_RUNS/MAFFT_ALL_RSVA.fasta)
MAFFT_ALL_RSVB_COUNT=$(grep -c "^>" $ALL_RUNS/MAFFT_ALL_RSVB.fasta)

step_complete "6" "Success!! Aligned RSVA: $MAFFT_ALL_RSVA_COUNT sequences!! Aligned RSVB: $MAFFT_ALL_RSVB_COUNT sequences!!"

########################################
# Step 7: Pan-run Tree Construction with iqtree
########################################
echo "Step 7: Pan-run Tree Construction with iqtree..."

# RSVA Tree
iqtree -s $ALL_RUNS/MAFFT_ALL_RSVA.fasta \
       -m MFP -T AUTO \
       -bb 1000 -alrt 1000 \
       -pre $ALL_RUNS/TREE_MAFFT_ALL_RSVA \
       2> $ALL_RUNS/trees_logs/iqtree_mafft_all_rsva.log

# RSVB Tree
iqtree -s $ALL_RUNS/MAFFT_ALL_RSVB.fasta \
       -m MFP -T AUTO \
       -bb 1000 -alrt 1000 \
       -pre $ALL_RUNS/TREE_MAFFT_ALL_RSVB \
       2> $ALL_RUNS/trees_logs/iqtree_mafft_all_rsvb.log

step_complete "7" "Success!!"

echo "========================================"
echo "Pan-run tree analysis successfully completed!"
echo "Output files:"
echo "- RSVA Alignment: $ALL_RUNS/MAFFT_ALL_RSVA.fasta"
echo "- RSVB Alignment: $ALL_RUNS/MAFFT_ALL_RSVB.fasta"
echo "- RSVA Tree: $ALL_RUNS/TREE_MAFFT_ALL_RSVA.treefile"
echo "- RSVB Tree: $ALL_RUNS/TREE_MAFFT_ALL_RSVB.treefile"
echo "========================================"
