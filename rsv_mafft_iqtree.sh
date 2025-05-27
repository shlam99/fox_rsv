#!/bin/bash
THREADS=$(nproc --all)

# Function to display step completion
step_complete() {
    echo ""
    echo "========================================"
    echo "Step $1 complete: $2"
    echo "========================================"
    echo ""
}

START_TIME=$(date +%s)

########################################
# Step 0: Setup and dependency checks
########################################
echo "Step 0: Setup and dependency checks..."

# Make directories
mkdir -p trees/{pooled_consensus,mafft_consensus,trees_logs,RSVA_tree,RSVB_tree,RSVAD_tree,RSVBD_tree}

# Check required tools
command -v mafft >/dev/null 2>&1 || { echo >&2 "Error: MAFFT not found"; exit 1; }
command -v iqtree >/dev/null 2>&1 || { echo >&2 "Error: IQ-TREE not found"; exit 1; }

step_complete "0" "Directory structure created and dependencies verified"

########################################
# Step 5: Pooling all consensus sequences by RSV type
########################################
echo "Step 5: Pooling all all consensus sequences by RSV type..."

# Process each RSV type individually
for type in RSVA RSVB RSVAD RSVBD; do
    if ls ${type}_* 1> /dev/null 2>&1; then
        # Combine and sort by name
        {
            cat ${type}_* 2>/dev/null | \
            awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
            sort -k1,1 | \
            tr '\t' '\n' > "trees/pooled_consensus/pooled_${type}.fasta"
        } || { echo "Error processing ${type} files"; continue; }
        
        COUNT=$(grep -c "^>" "trees/pooled_consensus/pooled_${type}.fasta" 2>/dev/null || echo 0)
        echo "Combined and sorted ${type} files (${COUNT} sequences)"
    else
        echo "No ${type} files found - skipping"
        continue
    fi
done

step_complete "5" "Files pooled in trees/pooled_consensus/"

########################################
# Step 6: MAFFT Alignment of pooled consensus of each type with mafft
########################################
echo "Step 6: MAFFT Alignment of pooled consensus of each type with mafft..."

for type in RSVA RSVB RSVAD RSVBD; do
    if [ -f "trees/pooled_consensus/pooled_${type}.fasta" ]; then
        echo "Aligning ${type}..."
        mafft --auto --thread $THREADS \
              --reorder \
              "trees/pooled_consensus/pooled_${type}.fasta" > \
              "trees/mafft_consensus/mafft_${type}.fasta" 2> \
              "trees/trees_logs/mafft_${type,,}.log"
        
        COUNT=$(grep -c "^>" "trees/mafft_consensus/mafft_${type}.fasta" 2>/dev/null || echo 0)
        echo "${type} aligned with ${COUNT} sequences"
    else
        echo "No ${type} pooled file found - skipping alignment"
    fi
done

MAFFT_RSVA_COUNT=$(grep -c "^>" trees/mafft_consensus/mafft_RSVA.fasta)
MAFFT_RSVB_COUNT=$(grep -c "^>" trees/mafft_consensus/mafft_RSVB.fasta)
MAFFT_RSVAD_COUNT=$(grep -c "^>" trees/mafft_consensus/mafft_RSVAD.fasta)
MAFFT_RSVBD_COUNT=$(grep -c "^>" trees/mafft_consensus/mafft_RSVBD.fasta)

step_complete "6" "Mafft Success!! Aligned RSVA: $MAFFT_RSVA_COUNT sequences!! Aligned RSVB: $MAFFT_RSVB_COUNT sequences!! Aligned RSVAD: $MAFFT_RSVAD_COUNT sequences!! Aligned RSVBD: $MAFFT_RSVBD_COUNT sequences!!"

echo ""
echo "========================================"
echo "Output files:"
echo "- RSVA Alignment: trees/mafft_consensus/mafft_RSVA.fasta"
echo "- RSVB Alignment: trees/mafft_consensus/mafft_RSVB.fasta"
echo "- RSVAD Alignment: trees/mafft_consensus/mafft_RSVAD.fasta"
echo "- RSVBD Alignment: trees/mafft_consensus/mafft_RSVBD.fasta"
echo "========================================"


########################################
# Step 7: Phylogenetic Tree Construction using iqtree
########################################
echo "Step 7: Phylogenetic Tree Construction using iqtree..."

for type in RSVA RSVB RSVAD RSVBD; do
    if [ -f "trees/mafft_consensus/mafft_${type}.fasta" ]; then
        echo "Building ${type} tree..."
        mkdir -p "trees/${type}_tree"
        iqtree -s "trees/mafft_consensus/mafft_${type}.fasta" \
               -m MFP -T $THREADS \
               -bb 1000 -alrt 1000 \
               -pre "trees/${type}_tree/${type}_tree" \
               2> "trees/trees_logs/iqtree_${type,,}.log"
        
        echo "${type} tree built in trees/${type}_tree/"
    else
        echo "No ${type} alignment found - skipping tree construction"
    fi
done

step_complete "7" "Tree Success!!"

echo ""
echo "========================================"
echo "Output files:"
echo "- RSVA Tree: trees/RSVA_tree/TREE_mafft_RSVA.treefile"
echo "- RSVB Tree: trees/RSVB_tree/TREE_mafft_RSVB.treefile"
echo "- RSVAD Tree: trees/RSVAD_tree/TREE_mafft_RSVAD.treefile"
echo "- RSVBD Tree: trees/RSVBD_tree/TREE_mafft_RSVBD.treefile"
echo "========================================"

END_TIME=$(date +%s)
echo "Total pipeline runtime: $((END_TIME - START_TIME)) seconds"
