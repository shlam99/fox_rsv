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

# Check for input files and convert if needed
HAS_BAM=$(ls *.bam 2>/dev/null | wc -l)
HAS_FASTQ=$(ls *.fastq *.fastq.gz 2>/dev/null | wc -l)

if [ "$HAS_BAM" -gt 0 ]; then
    echo "Found BAM files, converting directly to compressed FASTQ..."
    for i in {1..24}; do
        BARCODE_PADDED=$(printf "%02d" "$i")
        BAM_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.bam"
        FASTQ_OUT="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
        
        if [ -f "$BAM_IN" ]; then
            echo "Converting ${BAM_IN} directly to ${FASTQ_OUT}..."
            samtools fastq "$BAM_IN" | gzip > "$FASTQ_OUT"
        fi
    done
elif [ "$HAS_FASTQ" -gt 0 ]; then
    echo "Found FASTQ files, compressing if needed..."
    # Compress any uncompressed FASTQ files
    for f in *.fastq; do
        if [ -f "$f" ]; then
            echo "Compressing ${f}..."
            gzip "$f"
        fi
    done
else
    echo "Error: No input files found (expected either .bam or .fastq/.fastq.gz files)"
    exit 1
fi

step_complete "0" "Setup and dependency checks"


########################################
# Step 1: Quality filtering with filtlong
########################################
echo "Step 1: Performing quality filtering with filtlong..."

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    FASTQ_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    
    if [ -f "$FASTQ_IN" ]; then
        echo "Processing ${FASTQ_IN}..."
        filtlong --min_length $MIN_LENGTH \
                 --keep_percent $KEEP_PERCENT \
                 --target_bases $TARGET_BASES \
                 "$FASTQ_IN" 2> "qc_reads/qc_logs/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtlong.log" | \
        gzip > "$FASTQ_FILTERED"
    else
        echo "Warning: Input FASTQ file ${FASTQ_IN} not found"
    fi
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

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    IRMA_OUTDIR="irma_results/barcode${BARCODE_PADDED}"
    
    if [ -f "$FASTQ_FILTERED" ]; then
        echo "Processing ${FASTQ_FILTERED} with IRMA..."
        IRMA RSV_ont \
            "$FASTQ_FILTERED" \
            "$IRMA_OUTDIR"
    else
        echo "Warning: Filtered FASTQ file ${FASTQ_FILTERED} not found"
    fi
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

    # Add to appropriate consensus file with sample ID
    SAMPLE_ID="${SAMPLE_IDS[$BARCODE_ID]}"
    if [ -z "$SAMPLE_ID" ]; then
        echo "WARNING: No sample ID mapped for ${BARCODE_ID}" >&2
        SAMPLE_ID="${BARCODE_ID}"
    fi
    
    echo "Adding ${SAMPLE_ID} to ${RSV_TYPE} consensus..."
    sed "s/^>.*/>${SAMPLE_ID}/" "$CONSENSUS_SOURCE" >> "irma_consensus/${RSV_TYPE}_consensus.fasta"
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
# Step 4: Create batch-labeled versions
########################################
echo "Step 4: Creating batch-labeled consensus files (${BATCH})..."

for type in RSVA RSVB RSVAD RSVBD; do
    INPUT="irma_consensus/${type}_consensus.fasta"
    OUTPUT="irma_consensus/${type}_${BATCH}.fasta"
    
    if [ -f "$INPUT" ] && [ -s "$INPUT" ]; then
        cp "$INPUT" "$OUTPUT"
        echo "Created batch-labeled: $OUTPUT"
    else
        echo "No ${type} consensus found to label"
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