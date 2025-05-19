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
mkdir -p qc_reads/qc_logs irma_results consensus_sequences

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
echo "Output files: irma_results/barcode${BARCODE_PADDED}/consensus.fasta"
echo "========================================"

########################################
# Step 3: Select and pool consensus sequences
########################################
echo "Step 3: Selecting and pooling consensus sequences by RSV type..."

> consensus_sequences/RSVA_consensus.fasta
> consensus_sequences/RSVB_consensus.fasta
> consensus_sequences/RSVAD_consensus.fasta
> consensus_sequences/RSVBC_consensus.fasta
> consensus_sequences/RSVBD_consensus.fasta

for i in {1..22}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    IRMA_OUTDIR="irma_results/barcode${BARCODE_PADDED}"
    CONSENSUS="${IRMA_OUTDIR}/consensus.fasta"
    
    echo "Processing ${CONSENSUS}..."
    
    # Determine which RSV type is present in the IRMA output
    RSV_TYPE=$(grep "^>" "$CONSENSUS" | cut -d'|' -f2 | head -n1)
    
    # Add barcode information to the header and append to appropriate file
    case $RSV_TYPE in
        "RSVA")
            sed "s/^>/>${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}|/" "$CONSENSUS" >> consensus_sequences/RSVA_consensus.fasta
            ;;
        "RSVB")
            sed "s/^>/>${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}|/" "$CONSENSUS" >> consensus_sequences/RSVB_consensus.fasta
            ;;
        "RSVAD")
            sed "s/^>/>${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}|/" "$CONSENSUS" >> consensus_sequences/RSVAD_consensus.fasta
            ;;
        "RSVBC")
            sed "s/^>/>${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}|/" "$CONSENSUS" >> consensus_sequences/RSVBC_consensus.fasta
            ;;
        "RSVBD")
            sed "s/^>/>${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}|/" "$CONSENSUS" >> consensus_sequences/RSVBD_consensus.fasta
            ;;
        *)
            echo "WARNING: Unknown RSV type in ${CONSENSUS}" >&2
            ;;
    esac
done

# Count the number of sequences in each file
RSVA_COUNT=$(grep -c "^>" consensus_sequences/RSVA_consensus.fasta)
RSVB_COUNT=$(grep -c "^>" consensus_sequences/RSVB_consensus.fasta)
RSVAD_COUNT=$(grep -c "^>" consensus_sequences/RSVAD_consensus.fasta)
RSVBC_COUNT=$(grep -c "^>" consensus_sequences/RSVBC_consensus.fasta)
RSVBD_COUNT=$(grep -c "^>" consensus_sequences/RSVBD_consensus.fasta)

step_complete "3" "Success!!: Selected RSVs (A:$RSVA_COUNT B:$RSVB_COUNT AD:$RSVAD_COUNT BC:$RSVBC_COUNT BD:$RSVBD_COUNT)"

echo ""
echo "========================================"
echo "Output files:"
echo "- RSVA consensus: consensus_sequences/RSVA_consensus.fasta"
echo "- RSVB consensus: consensus_sequences/RSVB_consensus.fasta"
echo "- RSVAD consensus: consensus_sequences/RSVAD_consensus.fasta"
echo "- RSVBC consensus: consensus_sequences/RSVBC_consensus.fasta"
echo "- RSVBD consensus: consensus_sequences/RSVBD_consensus.fasta"
echo "========================================"