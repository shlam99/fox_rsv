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
# Step 0: Setup and dependency checks
########################################
echo "Starting pipeline..."
echo "Step 0: Setting up directory structure, checking input files, and checking dependencies..."

# Make directories
mkdir -p qc_reads/qc_logs irma_results irma_consensus

# Check for required tools
command -v samtools >/dev/null 2>&1 || { echo >&2 "Error: samtools not found."; exit 1; }
command -v filtlong >/dev/null 2>&1 || { echo >&2 "Error: filtlong not found."; exit 1; }
command -v IRMA >/dev/null 2>&1 || { echo >&2 "Error: IRMA not found."; exit 1; }

# Convert BAM to FASTQ (if needed)
if ls *.bam >/dev/null 2>&1; then
    echo "Found BAM files, converting to FASTQ in parallel..."
    for i in {1..24}; do
        BARCODE_PADDED=$(printf "%02d" "$i")
        BAM_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.bam"
        FASTQ_OUT="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
        
        if [ -f "$BAM_IN" ]; then
            samtools fastq "$BAM_IN" | gzip > "$FASTQ_OUT" &
        fi
        
        # Limit concurrent jobs to $THREADS (from config.sh)
        if [[ $(jobs -r -p | wc -l) -ge $THREADS ]]; then
            wait -n
        fi
    done
    wait  # Ensure all jobs finish
elif ls *.fastq >/dev/null 2>&1; then
    echo "Found FASTQ files, compressing in parallel..."
    find . -maxdepth 1 -name "*.fastq" -print0 | xargs -0 -P $THREADS -I {} sh -c 'echo "Compressing {}..."; gzip {}'
else
    echo "Error: No input files found (.bam or .fastq)"
fi

step_complete "0" "Setup and dependency checks"

########################################
# Step 1: Quality filtering with filtlong (parallelized)
########################################
echo "Step 1: Running filtlong in parallel (xargs)..."
seq 1 24 | xargs -P $THREADS -I {} bash -c '
    BARCODE_PADDED=$(printf "%02d" "$1")
    FASTQ_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    
    if [ -f "$FASTQ_IN" ]; then
        echo "Processing $FASTQ_IN..."
        filtlong --min_length $MIN_LENGTH \
                 --keep_percent $KEEP_PERCENT \
                 --target_bases $TARGET_BASES \
                 "$FASTQ_IN" 2> "qc_reads/qc_logs/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtlong.log" | \
        gzip > "$FASTQ_FILTERED"
    else
        echo "Warning: $FASTQ_IN not found" >&2
    fi
' _ {}

step_complete "1" "Quality filtering complete"

echo ""
echo "========================================"
echo "Output files: qc_reads/{SAMPLE_PREFIX}_barcode{BARCODE_PADDED}.filtered.fastq.gz"
echo "========================================"

########################################
# Step 2: Run IRMA in parallel (background jobs)
########################################
echo "Step 2: Running IRMA with $THREADS parallel jobs..."
for i in {1..24}; do
    (
        BARCODE_PADDED=$(printf "%02d" "$i")
        FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
        IRMA_OUTDIR="irma_results/barcode${BARCODE_PADDED}"
        
        if [ -f "$FASTQ_FILTERED" ]; then
            echo "Running IRMA on $FASTQ_FILTERED..."
            IRMA RSV_ont "$FASTQ_FILTERED" "$IRMA_OUTDIR"
        fi
    ) &
    
    # Limit concurrent IRMA jobs (RAM/CPU heavy)
    if [[ $(jobs -r -p | wc -l) -ge $THREADS ]]; then
        wait -n
    fi
done
wait  # Wait for all IRMA jobs

step_complete "2" "IRMA analysis complete"

########################################
# Step 3: Select and pool consensus sequences
########################################
echo "Step 3: Selecting and pooling consensus sequences by RSV type..."

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

# Count sequences
RSVA_COUNT=0; [ -f "irma_consensus/RSVA_consensus.fasta" ] && RSVA_COUNT=$(grep -c "^>" irma_consensus/RSVA_consensus.fasta)
RSVB_COUNT=0; [ -f "irma_consensus/RSVB_consensus.fasta" ] && RSVB_COUNT=$(grep -c "^>" irma_consensus/RSVB_consensus.fasta)
RSVAD_COUNT=0; [ -f "irma_consensus/RSVAD_consensus.fasta" ] && RSVAD_COUNT=$(grep -c "^>" irma_consensus/RSVAD_consensus.fasta)
RSVBD_COUNT=0; [ -f "irma_consensus/RSVBD_consensus.fasta" ] && RSVBD_COUNT=$(grep -c "^>" irma_consensus/RSVBD_consensus.fasta)

step_complete "3" "Success!!: Pooled RSVs (A:$RSVA_COUNT B:$RSVB_COUNT AD:$RSVAD_COUNT BD:$RSVBD_COUNT)"

echo ""
echo "========================================"
echo "Step 3 output files:"
echo "- RSVA: irma_consensus/RSVA_consensus.fasta"
echo "- RSVB: irma_consensus/RSVB_consensus.fasta"
echo "- RSVAD: irma_consensus/RSVAD_consensus.fasta"
echo "- RSVBD: irma_consensus/RSVBD_consensus.fasta"
echo "========================================"

########################################
# Step 4: Create batch-labeled versions with sample IDs
########################################
echo "Step 4: Creating batch-labeled consensus files (${BATCH})..."

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

# Count sequences in batch-labeled files
BATCH_RSVA_COUNT=$(grep -c "^>" irma_consensus/RSVA_${BATCH}.fasta 2>/dev/null || echo 0)
BATCH_RSVB_COUNT=$(grep -c "^>" irma_consensus/RSVB_${BATCH}.fasta 2>/dev/null || echo 0)
BATCH_RSVAD_COUNT=$(grep -c "^>" irma_consensus/RSVAD_${BATCH}.fasta 2>/dev/null || echo 0)
BATCH_RSVBD_COUNT=$(grep -c "^>" irma_consensus/RSVBD_${BATCH}.fasta 2>/dev/null || echo 0)

step_complete "4" "Success!!: Batch-labeled RSVs (A:$BATCH_RSVA_COUNT B:$BATCH_RSVB_COUNT AD:$BATCH_RSVAD_COUNT BD:$BATCH_RSVBD_COUNT)"

echo ""
echo "========================================"
echo "Step 4 output files:"
echo "- RSVA: irma_consensus/RSVA_${BATCH}.fasta"
echo "- RSVB: irma_consensus/RSVB_${BATCH}.fasta"
echo "- RSVAD: irma_consensus/RSVAD_${BATCH}.fasta"
echo "- RSVBD: irma_consensus/RSVBD_${BATCH}.fasta"
echo "========================================"

END_TIME=$(date +%s)
echo "Total pipeline runtime: $((END_TIME - START_TIME)) seconds"

