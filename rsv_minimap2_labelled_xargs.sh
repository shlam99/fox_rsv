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

# Make directories (including pooled_consensus)
mkdir -p qc_reads/qc_logs aligned_bams consensus_sequences pooled_consensus

# Check for required tools
command -v samtools >/dev/null 2>&1 || { echo >&2 "Error: samtools not found."; exit 1; }
command -v filtlong >/dev/null 2>&1 || { echo >&2 "Error: filtlong not found."; exit 1; }

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
        
        # Limit concurrent jobs to $THREADS
        if [[ $(jobs -r -p | wc -l) -ge $THREADS ]]; then
            wait -n
        fi
    done
    wait
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
echo "Output files: qc_reads/${SAMPLE_PREFIX}_barcode01-24.filtered.fastq.gz"
echo "========================================"

########################################
# Step 2: Alignment with minimap2
########################################
echo "Step 2: Aligning reads with minimap2 in parallel (xargs)..."

seq 1 24 | xargs -P $THREADS -I {} bash -c '
    BARCODE_PADDED=$(printf "%02d" "$1")
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    BAM_PREFIX="aligned_bams/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}"
    
    if [ -f "$FASTQ_FILTERED" ]; then
        echo "Aligning $FASTQ_FILTERED..."
        minimap2 -ax map-ont \
                 -t 1 \
                 --MD \
                 -Y \
                 "$REFERENCE_GENOME" \
                 "$FASTQ_FILTERED" | \
        samtools sort -@ 1 -o "${BAM_PREFIX}.aligned.bam"
        samtools index -@ 1 "${BAM_PREFIX}.aligned.bam"
        
        # Generate alignment stats
        samtools flagstat "${BAM_PREFIX}.aligned.bam" > "${BAM_PREFIX}.flagstat"
        samtools stats "${BAM_PREFIX}.aligned.bam" > "${BAM_PREFIX}.stats"
    else
        echo "Warning: $FASTQ_FILTERED not found" >&2
    fi
' _ {}

step_complete "2" "Alignment complete"

echo ""
echo "========================================"
echo "Output files: aligned_bams/${SAMPLE_PREFIX}_barcode01-24.aligned.bam"
echo "========================================"

########################################
# Step 3: Consensus generation with samtools
########################################
echo "Step 3: Generating consensus sequences in parallel (xargs)..."
seq 1 24 | xargs -P $THREADS -I {} bash -c '
    BARCODE_PADDED=$(printf "%02d" "$1")
    BAM="aligned_bams/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.aligned.bam"
    CONSENSUS="consensus_sequences/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fasta"
    
    if [ -f "$BAM" ]; then
        echo "Processing $BAM..."
        
        # Verify BAM is aligned to this reference
        if ! samtools view -H "$BAM" | grep -q "$(head -n1 "$REFERENCE_GENOME" | tr -d ">")"; then
            echo "ERROR: BAM not aligned to specified reference" >&2
            exit 1
        fi
        
        samtools consensus -a --show-ins no "$BAM" -o "$CONSENSUS"
        
        if [ ! -s "$CONSENSUS" ]; then
            echo "ERROR: Empty consensus for barcode ${BARCODE_PADDED}" >&2
            exit 1
        fi
    else
        echo "Warning: $BAM not found" >&2
    fi
' _ {}

step_complete "3" "Consensus generation complete"

echo ""
echo "========================================"
echo "Output files: consensus_sequences/${SAMPLE_PREFIX}_barcode01-24.fasta"
echo "========================================"

########################################
# Step 4: Select RSVA OR RSVB consensus per sample (with barcode+sample ID)
########################################
echo "Step 4: Select highest-quality RSVA/RSVB consensus per sample (with barcode+sample ID)..."

# Initialize output files
> "pooled_consensus/RSVA_consensus_${BATCH}.fasta"
> "pooled_consensus/RSVB_consensus_${BATCH}.fasta"

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BARCODE_ID="barcode${BARCODE_PADDED}"
    INPUT_FASTA="consensus_sequences/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fasta"
    
    echo "Processing ${INPUT_FASTA}..."
    
    # Get sample ID
    SAMPLE_ID=$(get_sample_id "$BARCODE_ID")
    
    # Count non-N bases for each virus in the sample
    RSVA_COUNT=$(awk '/^>RSVA/{getline; seq=""; while ((getline line)>0) {if(line~/^>/){break}; seq=seq line}; gsub(/[^ACGT]/, "", seq); print length(seq)}' "$INPUT_FASTA")
    RSVB_COUNT=$(awk '/^>RSVB/{getline; seq=""; while ((getline line)>0) {if(line~/^>/){break}; seq=seq line}; gsub(/[^ACGT]/, "", seq); print length(seq)}' "$INPUT_FASTA")
    
    # Select the version with more ACGT bases (fewer Ns)
    if [ "$RSVA_COUNT" -ge "$RSVB_COUNT" ]; then
        # Use RSVA version with barcode and sample ID
        awk -v batch="$BATCH" -v bar="$BARCODE_PADDED" -v sample="$SAMPLE_ID" '
            /^>RSVA/ { 
                print $0 "|" batch "_barcode" bar "|" sample
                while ((getline line) > 0) {
                    if (line ~ /^>/) exit
                    print line
                }
            }
        ' "$INPUT_FASTA" >> "pooled_consensus/RSVA_consensus_${BATCH}.fasta"
    else
        # Use RSVB version with barcode and sample ID
        awk -v batch="$BATCH" -v bar="$BARCODE_PADDED" -v sample="$SAMPLE_ID" '
            /^>RSVB/ { 
                print $0 "|" batch "_barcode" bar "|" sample
                while ((getline line) > 0) {
                    if (line ~ /^>/) exit
                    print line
                }
            }
        ' "$INPUT_FASTA" >> "pooled_consensus/RSVB_consensus_${BATCH}.fasta"
    fi
done

RSVA_COUNT=$(grep -c "^>RSVA" "pooled_consensus/RSVA_consensus_${BATCH}.fasta")
RSVB_COUNT=$(grep -c "^>RSVB" "pooled_consensus/RSVB_consensus_${BATCH}.fasta")

step_complete "4" "Success!!: Selected RSVA ($RSVA_COUNT samples) and RSVB ($RSVB_COUNT samples) with barcode+sample ID"

echo ""
echo "========================================"
echo "Step 4 output files:"
echo "- RSVA: pooled_consensus/RSVA_consensus_${BATCH}.fasta"
echo "- RSVB: pooled_consensus/RSVB_consensus_${BATCH}.fasta"
echo "========================================"

########################################
# Step 5: Create simplified consensus with sample ID only
########################################
echo "Step 5: Creating simplified consensus files with sample ID only..."

# Initialize output files
> "pooled_consensus/RSVA_${BATCH}.fasta"
> "pooled_consensus/RSVB_${BATCH}.fasta"

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BARCODE_ID="barcode${BARCODE_PADDED}"
    INPUT_FASTA="consensus_sequences/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fasta"
    
    echo "Processing ${INPUT_FASTA} for simplified consensus..."
    
    # Get sample ID
    SAMPLE_ID=$(get_sample_id "$BARCODE_ID")
    
    # Count non-N bases for each virus in the sample
    RSVA_COUNT=$(awk '/^>RSVA/{getline; seq=""; while ((getline line)>0) {if(line~/^>/){break}; seq=seq line}; gsub(/[^ACGT]/, "", seq); print length(seq)}' "$INPUT_FASTA")
    RSVB_COUNT=$(awk '/^>RSVB/{getline; seq=""; while ((getline line)>0) {if(line~/^>/){break}; seq=seq line}; gsub(/[^ACGT]/, "", seq); print length(seq)}' "$INPUT_FASTA")
    
    # Select the version with more ACGT bases (fewer Ns)
    if [ "$RSVA_COUNT" -ge "$RSVB_COUNT" ]; then
        # Use RSVA version with sample ID only
        awk -v sample="$SAMPLE_ID" '
            /^>RSVA/ { 
                print ">" sample
                while ((getline line) > 0) {
                    if (line ~ /^>/) exit
                    print line
                }
            }
        ' "$INPUT_FASTA" >> "pooled_consensus/RSVA_${BATCH}.fasta"
    else
        # Use RSVB version with sample ID only
        awk -v sample="$SAMPLE_ID" '
            /^>RSVB/ { 
                print ">" sample
                while ((getline line) > 0) {
                    if (line ~ /^>/) exit
                    print line
                }
            }
        ' "$INPUT_FASTA" >> "pooled_consensus/RSVB_${BATCH}.fasta"
    fi
done

RSVA_COUNT=$(grep -c "^>" "pooled_consensus/RSVA_${BATCH}.fasta")
RSVB_COUNT=$(grep -c "^>" "pooled_consensus/RSVB_${BATCH}.fasta")

step_complete "5" "Success!!: Created simplified RSVA ($RSVA_COUNT samples) and RSVB ($RSVB_COUNT samples) with sample ID only"

echo ""
echo "========================================"
echo "Step 5 output files:"
echo "- RSVA: pooled_consensus/RSVA_${BATCH}.fasta"
echo "- RSVB: pooled_consensus/RSVB_${BATCH}.fasta"
echo "========================================"

END_TIME=$(date +%s)
echo "Total pipeline runtime: $((END_TIME - START_TIME)) seconds"
