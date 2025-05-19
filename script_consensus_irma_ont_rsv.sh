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
echo "Starting pipeline..."
echo "Step 0: Setting up directory structure and checking dependencies..."

# Make directories
mkdir -p qc_reads/qc_logs irma_results irma_consensus

# Check for required tools
command -v filtlong >/dev/null 2>&1 || { echo >&2 "Error: filtlong not found."; exit 1; }
command -v IRMA >/dev/null 2>&1 || { echo >&2 "Error: IRMA not found."; exit 1; }

# Compress FASTQ files if needed
if ls *.fastq >/dev/null 2>&1; then
    echo "Compressing FASTQ files..."
    gzip *.fastq
fi

step_complete "0" "Setup and dependency checks"

########################################
# Step 1: Quality filtering with filtlong
########################################
echo "Step 1: Performing quality filtering with filtlong..."

for i in {1..22}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    FASTQ_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    
    echo "Processing ${FASTQ_IN}..."
    filtlong --min_length $MIN_LENGTH \
             --keep_percent $KEEP_PERCENT \
             --target_bases $TARGET_BASES \
             "$FASTQ_IN" 2> "qc_reads/qc_logs/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtlong.log" | \
    gzip > "$FASTQ_FILTERED"
done

step_complete "1" "Success!!"

echo ""
echo "========================================"
echo "Output files: qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
echo "========================================"

########################################
# Step 2: Run IRMA for each barcode
########################################
echo "Step 2: Running IRMA for each barcode..."

for i in {1..22}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    IRMA_OUTDIR="irma_results/barcode${BARCODE_PADDED}"
    
    echo "Processing ${FASTQ_FILTERED} with IRMA..."
    
    IRMA RSV_ont \
        "$FASTQ_FILTERED" \
        "$IRMA_OUTDIR" \

done

step_complete "2" "Success!!"

echo ""
echo "========================================"
echo "Output files:"
echo "- Consensus sequence: irma_results/amended_sequences/barcode${BARCODE_PADDED}.fa"
echo "- Reference used: irma_results/RSV_*.fasta"
echo "========================================"

########################################
# Step 3: Select and pool consensus sequences
########################################
echo "Step 3: Selecting and pooling consensus sequences by RSV type..."

# Initialize output files
> irma_consensus/RSVA_consensus.fasta
> irma_consensus/RSVB_consensus.fasta
> irma_consensus/RSVAD_consensus.fasta
> irma_consensus/RSVBD_consensus.fasta

for i in {1..22}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    IRMA_OUTDIR="irma_results/barcode${BARCODE_PADDED}"
    CONSENSUS_SOURCE="${IRMA_OUTDIR}/amended_consensus/barcode${BARCODE_PADDED}.fa"
    TYPE_FILE="${IRMA_OUTDIR}/RSV_"*".fasta"  # This will match RSV_A.fasta, RSV_B.fasta, etc.
    
    echo "Processing barcode ${BARCODE_PADDED}..."

    # Check if consensus file exists
    if [ ! -f "$CONSENSUS_SOURCE" ]; then
        echo "WARNING: Consensus file not found for barcode ${BARCODE_PADDED}" >&2
        continue
    fi

    # Determine RSV type from the parent directory files
    RSV_TYPE=""
    for file in $TYPE_FILE; do
        case $(basename "$file") in
            "RSV_A.fasta") RSV_TYPE="RSVA" ;;
            "RSV_B.fasta") RSV_TYPE="RSVB" ;;
            "RSV_AD.fasta") RSV_TYPE="RSVAD" ;;
            "RSV_BD.fasta") RSV_TYPE="RSVBD" ;;
        esac
        break  # We only need to check one file
    done

    if [ -z "$RSV_TYPE" ]; then
        echo "WARNING: Could not determine RSV type for barcode ${BARCODE_PADDED}" >&2
        continue
    fi

    # Add to appropriate consensus file
    echo "Adding barcode ${BARCODE_PADDED} to ${RSV_TYPE} consensus..."
    sed "s/^>/>${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}|/" "$CONSENSUS_SOURCE" >> "irma_consensus/${RSV_TYPE}_consensus.fasta"
done

# Count the number of sequences in each file
RSVA_COUNT=$(grep -c "^>" irma_consensus/RSVA_consensus.fasta || echo 0)
RSVB_COUNT=$(grep -c "^>" irma_consensus/RSVB_consensus.fasta || echo 0)
RSVAD_COUNT=$(grep -c "^>" irma_consensus/RSVAD_consensus.fasta || echo 0)
RSVBD_COUNT=$(grep -c "^>" irma_consensus/RSVBD_consensus.fasta || echo 0)

step_complete "3" "Success!!: Selected RSVs (A:$RSVA_COUNT B:$RSVB_COUNT AD:$RSVAD_COUNT BD:$RSVBD_COUNT)"

echo ""
echo "========================================"
echo "Output files:"
echo "- RSVA consensus: irma_consensus/RSVA_consensus.fasta"
echo "- RSVB consensus: irma_consensus/RSVB_consensus.fasta"
echo "- RSVAD consensus: irma_consensus/RSVAD_consensus.fasta"
echo "- RSVBD consensus: irma_consensus/RSVBD_consensus.fasta"
echo "========================================"
