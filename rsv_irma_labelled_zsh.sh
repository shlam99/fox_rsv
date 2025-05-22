#!/usr/bin/env bash
source config.sh

# Function to display step completion
step_complete() {
    printf "\n========================================\n"
    printf "Step %s complete: %s\n" "$1" "$2"
    printf "========================================\n\n"
}

START_TIME=$(date +%s)

########################################
# Step 0: Setup and dependency checks
########################################
echo "Starting pipeline..."
echo "Step 0: Setting up directories and checking dependencies..."

mkdir -p qc_reads/qc_logs irma_results irma_consensus

# Check tools (compatible syntax)
for cmd in samtools filtlong IRMA; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo >&2 "Error: $cmd not found."
        exit 1
    fi
done

# Convert BAM to FASTQ (parallel)
if ls *.bam >/dev/null 2>&1; then
    echo "Converting BAM to FASTQ..."
    for i in $(seq 1 24); do
        BARCODE_PADDED=$(printf "%02d" "$i")
        BAM_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.bam"
        FASTQ_OUT="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
        
        if [ -f "$BAM_IN" ]; then
            samtools fastq "$BAM_IN" | gzip > "$FASTQ_OUT" &
        fi
        
        # Limit jobs to THREADS
        if [[ $(jobs -r -p | wc -l) -ge $THREADS ]]; then
            wait -n || true
        fi
    done
    wait
elif ls *.fastq >/dev/null 2>&1; then
    echo "Compressing FASTQ files..."
    find . -maxdepth 1 -name "*.fastq" -print0 | \
        xargs -0 -P $THREADS -I {} sh -c 'gzip -f "{}"'
else
    echo "Error: No input files found (.bam or .fastq)"
fi

step_complete "0" "Setup complete"

########################################
# Step 1: Quality filtering (parallel xargs)
########################################
echo "Running filtlong with $THREADS threads..."
seq 1 24 | xargs -P $THREADS -I {} bash -c '
    BARCODE_PADDED=$(printf "%02d" "{}")
    FASTQ_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    
    if [ -f "$FASTQ_IN" ]; then
        filtlong --min_length $MIN_LENGTH \
                 --keep_percent $KEEP_PERCENT \
                 --target_bases $TARGET_BASES \
                 "$FASTQ_IN" 2> "qc_reads/qc_logs/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtlong.log" | \
        gzip > "$FASTQ_FILTERED"
    fi
'

step_complete "1" "Quality filtering done"

########################################
# Step 2: IRMA (background jobs)
########################################
echo "Running IRMA with $THREADS parallel jobs..."
for i in $(seq 1 24); do
    (
        BARCODE_PADDED=$(printf "%02d" "$i")
        FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
        IRMA_OUTDIR="irma_results/barcode${BARCODE_PADDED}"
        
        if [ -f "$FASTQ_FILTERED" ]; then
            IRMA RSV_ont "$FASTQ_FILTERED" "$IRMA_OUTDIR"
        fi
    ) &
    
    # Job control
    if [[ $(jobs -r -p | wc -l) -ge $THREADS ]]; then
        wait -n || true
    fi
done
wait

step_complete "2" "IRMA analysis done"

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
    TYPE_FILE=$(ls ${IRMA_OUTDIR}/RSV_*.fasta 2>/dev/null | head -1)
    
    if [ ! -f "$CONSENSUS_SOURCE" ]; then
        echo "WARNING: Consensus file not found for ${BARCODE_ID}" >&2
        continue
    fi

    # Determine RSV type
    RSV_TYPE=""
    case $(basename "$TYPE_FILE") in
        "RSV_A.fasta") RSV_TYPE="RSVA" ;;
        "RSV_B.fasta") RSV_TYPE="RSVB" ;;
        "RSV_AD.fasta") RSV_TYPE="RSVAD" ;;
        "RSV_BD.fasta") RSV_TYPE="RSVBD" ;;
        *) echo "WARNING: Could not determine RSV type for ${BARCODE_ID}" >&2 ;;
    esac

    if [ -n "$RSV_TYPE" ]; then
        SAMPLE_ID="${SAMPLE_IDS[$BARCODE_ID]}"
        if [ -z "$SAMPLE_ID" ]; then
            SAMPLE_ID="${BARCODE_ID}"
        fi
        echo "Adding ${SAMPLE_ID} to ${RSV_TYPE} consensus..."
        sed "s/^>.*/>${SAMPLE_ID}/" "$CONSENSUS_SOURCE" >> "irma_consensus/${RSV_TYPE}_consensus.fasta"
    fi
done

# Remove empty consensus files
for type in RSVA RSVB RSVAD RSVBD; do
    FILE="irma_consensus/${type}_consensus.fasta"
    if [ -f "$FILE" ] && [ ! -s "$FILE" ]; then
        rm "$FILE"
    fi
done

# Count sequences
RSVA_COUNT=$(grep -c "^>" irma_consensus/RSVA_consensus.fasta 2>/dev/null || echo 0)
RSVB_COUNT=$(grep -c "^>" irma_consensus/RSVB_consensus.fasta 2>/dev/null || echo 0)
RSVAD_COUNT=$(grep -c "^>" irma_consensus/RSVAD_consensus.fasta 2>/dev/null || echo 0)
RSVBD_COUNT=$(grep -c "^>" irma_consensus/RSVBD_consensus.fasta 2>/dev/null || echo 0)

step_complete "3" "Pooled RSVs (A:$RSVA_COUNT B:$RSVB_COUNT AD:$RSVAD_COUNT BD:$RSVBD_COUNT)"

########################################
# Step 4: Create batch-labeled versions
########################################
echo "Step 4: Creating batch-labeled consensus files (${BATCH})..."

for type in RSVA RSVB RSVAD RSVBD; do
    INPUT="irma_consensus/${type}_consensus.fasta"
    OUTPUT="irma_consensus/${type}_${BATCH}.fasta"
    
    if [ -f "$INPUT" ] && [ -s "$INPUT" ]; then
        cp "$INPUT" "$OUTPUT"
        echo "Created batch-labeled: $OUTPUT"
    fi
done

step_complete "4" "Batch labeling complete"

END_TIME=$(date +%s)
echo "Total pipeline runtime: $((END_TIME - START_TIME)) seconds"
