#!/bin/bash
source config.sh

# Export all variables needed for parallel processing
export SAMPLE_PREFIX BATCH THREADS MIN_LENGTH TARGET_BASES KEEP_PERCENT

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
# Step 3: Pool typed consensus sequences into BatchID/BarcodeXX/SampleID labelled consensus
########################################
echo "Step 3: Pool typed consensus sequences into BatchID/BarcodeXX/SampleID labelled consensus for ${BATCH}..."

# Make directory
mkdir -p irma_consensus

# Initialize output files
> irma_consensus/RSVA_consensus.fasta
> irma_consensus/RSVB_consensus.fasta
> irma_consensus/RSVAD_consensus.fasta
> irma_consensus/RSVBD_consensus.fasta

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BARCODE_ID="barcode${BARCODE_PADDED}"
    IRMA_OUTDIR="irma_results/${BARCODE_ID}"
    CONSENSUS_SOURCE="${IRMA_OUTDIR}/amended_consensus/${BARCODE_ID}.fa"
    TYPE_FILE="${IRMA_OUTDIR}/RSV_"*".fasta"
    
    echo "Processing ${BARCODE_ID}..."

    # Check if consensus file exists
    if [ ! -f "$CONSENSUS_SOURCE" ]; then
        echo "WARNING: Consensus file not found for ${BARCODE_ID}" >&2
        continue
    fi

    # Determine RSV type
    RSV_TYPE=""
    for file in $TYPE_FILE; do
        case $(basename "$file") in
            "RSV_A.fasta") RSV_TYPE="RSVA" ;;
            "RSV_B.fasta") RSV_TYPE="RSVB" ;;
            "RSV_AD.fasta") RSV_TYPE="RSVAD" ;;
            "RSV_BD.fasta") RSV_TYPE="RSVBD" ;;
        esac
        break
    done

    if [ -z "$RSV_TYPE" ]; then
        echo "WARNING: Could not determine RSV type for ${BARCODE_ID}" >&2
        continue
    fi

    # Get sample ID for batch file 
    SAMPLE_ID=$(get_sample_id "$BARCODE_ID")

    # Add to pooled consensus with barcode prefix + sample ID (Step 3 format)
    sed "s/^>/>${BATCH}_barcode${BARCODE_PADDED}|${SAMPLE_ID}|/" "$CONSENSUS_SOURCE" >> "irma_consensus/${RSV_TYPE}_consensus_${BATCH}.fasta"
done

# Remove empty consensus files
for type in RSVA RSVB RSVAD RSVBD; do
    FILE="irma_consensus/${type}_consensus.fasta"
    if [ -f "$FILE" ] && [ ! -s "$FILE" ]; then
        echo "Removing empty file: $FILE"
        rm "$FILE"
    fi
done

# Count sequences in batch-labeled files
RSVA_COUNT=$(grep -c "^>" irma_consensus/RSVA_consensus_${BATCH}.fasta 2>/dev/null || echo 0)
RSVB_COUNT=$(grep -c "^>" irma_consensus/RSVB_consensus_${BATCH}.fasta 2>/dev/null || echo 0)
RSVAD_COUNT=$(grep -c "^>" irma_consensus/RSVAD_consensus_${BATCH}.fasta 2>/dev/null || echo 0)
RSVBD_COUNT=$(grep -c "^>" irma_consensus/RSVBD_consensus_${BATCH}.fasta 2>/dev/null || echo 0)

step_complete "3" "Pooled!! BatchID/BarcodeXX/SampleID RSVs (A:$RSVA_COUNT B:$RSVB_COUNT AD:$RSVAD_COUNT BD:$RSVBD_COUNT)"

echo ""
echo "========================================"
echo "Step 3 output files:"
echo "- RSVA: irma_consensus/RSVA_consensus_${BATCH}.fasta"
echo "- RSVB: irma_consensus/RSVB_consensus_${BATCH}.fasta"
echo "- RSVAD: irma_consensus/RSVAD_consensus_${BATCH}.fasta"
echo "- RSVBD: irma_consensus/RSVBD_consensus_${BATCH}.fasta"
echo "========================================"

########################################
# Step 4: Pool typed consensus sequences into SampleID labelled consensus
########################################
echo "Step 4: Pool typed consensus sequences into SampleID only labelled consensus for ${BATCH}..."

# Make directory
mkdir -p irma_consensus

# Initialize batch files
> irma_consensus/RSVA_${BATCH}.fasta
> irma_consensus/RSVB_${BATCH}.fasta
> irma_consensus/RSVAD_${BATCH}.fasta
> irma_consensus/RSVBD_${BATCH}.fasta

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BARCODE_ID="barcode${BARCODE_PADDED}"
    IRMA_OUTDIR="irma_results/${BARCODE_ID}"
    CONSENSUS_SOURCE="${IRMA_OUTDIR}/amended_consensus/${BARCODE_ID}.fa"
    TYPE_FILE="${IRMA_OUTDIR}/RSV_"*".fasta"
    
    # Check if consensus file exists
    if [ ! -f "$CONSENSUS_SOURCE" ]; then
        continue
    fi

    # Determine RSV type
    RSV_TYPE=""
    for file in $TYPE_FILE; do
        case $(basename "$file") in
            "RSV_A.fasta") RSV_TYPE="RSVA" ;;
            "RSV_B.fasta") RSV_TYPE="RSVB" ;;
            "RSV_AD.fasta") RSV_TYPE="RSVAD" ;;
            "RSV_BD.fasta") RSV_TYPE="RSVBD" ;;
        esac
        break
    done

    if [ -z "$RSV_TYPE" ]; then
        continue
    fi

    # Get sample ID for batch file 
    SAMPLE_ID=$(get_sample_id "$BARCODE_ID")
     
    # Add to pooled consensus with sample ID only (Step 4 format) 
    sed "s/^>.*/>${SAMPLE_ID}/" "$CONSENSUS_SOURCE" >> "irma_consensus/${RSV_TYPE}_${BATCH}.fasta" 
done 

# Remove empty consensus files
for type in RSVA RSVB RSVAD RSVBD; do
    FILE="irma_consensus/${type}_${BATCH}.fasta"
    if [ -f "$FILE" ] && [ ! -s "$FILE" ]; then
        echo "Removing empty file: $FILE"
        rm "$FILE"
    fi
done

# Count sequences in sampleID-labelled files
BATCH_RSVA_COUNT=$(grep -c "^>" irma_consensus/RSVA_${BATCH}.fasta 2>/dev/null || echo 0)
BATCH_RSVB_COUNT=$(grep -c "^>" irma_consensus/RSVB_${BATCH}.fasta 2>/dev/null || echo 0)
BATCH_RSVAD_COUNT=$(grep -c "^>" irma_consensus/RSVAD_${BATCH}.fasta 2>/dev/null || echo 0)
BATCH_RSVBD_COUNT=$(grep -c "^>" irma_consensus/RSVBD_${BATCH}.fasta 2>/dev/null || echo 0)

step_complete "4" "Pooled!! SampleID-labeled RSVs (A:$BATCH_RSVA_COUNT B:$BATCH_RSVB_COUNT AD:$BATCH_RSVAD_COUNT BD:$BATCH_RSVBD_COUNT)"

echo ""
echo "========================================"
echo "Step 4 output files:"
echo "- RSVA: irma_consensus/RSVA_${BATCH}.fasta"
echo "- RSVB: irma_consensus/RSVB_${BATCH}.fasta"
echo "- RSVAD: irma_consensus/RSVAD_${BATCH}.fasta"
echo "- RSVBD: irma_consensus/RSVBD_${BATCH}.fasta"
echo "========================================"

END_TIME=$(date +%s)
echo "Total pipeline runtime for ${BATCH}: $((END_TIME - START_TIME)) seconds"