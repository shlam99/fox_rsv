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
echo "Starting fox IAV pipeline on ${BATCH}..."
echo "Sample prefix is ${SAMPLE_PREFIX}..."
echo "Step 0: Setting up directory structure, checking input files, and checking dependencies..."

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

# Make directory
mkdir -p qc_reads/qc_logs

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
echo "Output files: qc_reads/${SAMPLE_PREFIX}_barcode01-24.filtered.fastq.gz"
echo "========================================"

########################################
# Step 2: Run IRMA in parallel (background jobs)
########################################
echo "Step 2: Running IRMA with $THREADS parallel jobs..."

# Make directory
mkdir -p irma_results

for i in {1..24}; do
    (
        BARCODE_PADDED=$(printf "%02d" "$i")
        FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
        IRMA_OUTDIR="irma_results/barcode${BARCODE_PADDED}"
        
        if [ -f "$FASTQ_FILTERED" ]; then
            echo "Running IRMA on $FASTQ_FILTERED..."
            IRMA FLU_ont "$FASTQ_FILTERED" "$IRMA_OUTDIR"
        fi
    ) &
    
    # Limit concurrent IRMA jobs (RAM/CPU heavy)
    if [[ $(jobs -r -p | wc -l) -ge $THREADS ]]; then
        wait -n
    fi
done
wait  # Wait for all IRMA jobs

step_complete "2" "IRMA analysis complete"

echo ""
echo "========================================"
echo "Output files: irma_results/barcode01-24/amended_consensus/barcode01-24.fa"
echo "========================================"

########################################
# Step 3: Pool consensus sequences into BatchID/BarcodeXX/SampleID labelled consensus
########################################
echo "Step 3: Pool consensus sequences into BatchID/BarcodeXX/SampleID labelled consensus for ${BATCH}..."

# Make directory
mkdir -p irma_consensus

# Initialize output file
> irma_consensus/IAV_consensus_${BATCH}.fasta

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BARCODE_ID="barcode${BARCODE_PADDED}"
    IRMA_OUTDIR="irma_results/${BARCODE_ID}"
    CONSENSUS_SOURCE="${IRMA_OUTDIR}/amended_consensus/${BARCODE_ID}.fa"
    
    echo "Processing ${BARCODE_ID}..."

    # Check if consensus file exists
    if [ ! -f "$CONSENSUS_SOURCE" ]; then
        echo "WARNING: Consensus file not found for ${BARCODE_ID}" >&2
        continue
    fi

    # Get sample ID for batch file 
    SAMPLE_ID=$(get_sample_id "$BARCODE_ID")

    # Add to pooled consensus with barcode prefix + sample ID (Step 3 format)
    sed "s/^>/>${BATCH}_barcode${BARCODE_PADDED}|${SAMPLE_ID}|/" "$CONSENSUS_SOURCE" >> "irma_consensus/IAV_consensus_${BATCH}.fasta"
done

# Remove empty consensus file if no sequences were added
if [ ! -s "irma_consensus/IAV_consensus_${BATCH}.fasta" ]; then
    echo "Removing empty file: irma_consensus/IAV_consensus_${BATCH}.fasta"
    rm "irma_consensus/IAV_consensus_${BATCH}.fasta"
fi

# Count sequences in batch-labeled file
IAV_COUNT=$(grep -c "^>" irma_consensus/IAV_consensus_${BATCH}.fasta 2>/dev/null || echo 0)

step_complete "3" "Pooled!! BatchID/BarcodeXX/SampleID IAV sequences ($IAV_COUNT total)"

echo ""
echo "========================================"
echo "Step 3 output file:"
echo "- IAV: irma_consensus/IAV_consensus_${BATCH}.fasta"
echo "========================================"

########################################
# Step 4: Pool consensus sequences into SampleID labelled consensus
########################################
echo "Step 4: Pool consensus sequences into SampleID only labelled consensus for ${BATCH}..."

# Make directory
mkdir -p irma_consensus

# Initialize batch file
> irma_consensus/IAV_${BATCH}.fasta

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BARCODE_ID="barcode${BARCODE_PADDED}"
    IRMA_OUTDIR="irma_results/${BARCODE_ID}"
    CONSENSUS_SOURCE="${IRMA_OUTDIR}/amended_consensus/${BARCODE_ID}.fa"
    
    # Check if consensus file exists
    if [ ! -f "$CONSENSUS_SOURCE" ]; then
        continue
    fi

    # Get sample ID for batch file 
    SAMPLE_ID=$(get_sample_id "$BARCODE_ID")
     
    # Add to pooled consensus with sample ID only (Step 4 format) 
    sed "s/^>.*/>${SAMPLE_ID}/" "$CONSENSUS_SOURCE" >> "irma_consensus/IAV_${BATCH}.fasta" 
done 

# Remove empty consensus file if no sequences were added
if [ ! -s "irma_consensus/IAV_${BATCH}.fasta" ]; then
    echo "Removing empty file: irma_consensus/IAV_${BATCH}.fasta"
    rm "irma_consensus/IAV_${BATCH}.fasta"
fi

# Count sequences in sampleID-labelled file
BATCH_IAV_COUNT=$(grep -c "^>" irma_consensus/IAV_${BATCH}.fasta 2>/dev/null || echo 0)

step_complete "4" "Pooled!! SampleID-labeled IAV sequences ($BATCH_IAV_COUNT total)"

echo ""
echo "========================================"
echo "Step 4 output file:"
echo "- IAV: irma_consensus/IAV_${BATCH}.fasta"
echo "========================================"

END_TIME=$(date +%s)
echo "Total pipeline runtime for ${BATCH}: $((END_TIME - START_TIME)) seconds"
