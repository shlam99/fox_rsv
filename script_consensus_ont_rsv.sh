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
mkdir -p qc_reads/qc_logs aligned_bams consensus_sequences 

# Check for required tools
command -v filtlong >/dev/null 2>&1 || { echo >&2 "Error: filtlong not found."; exit 1; }
command -v minimap2 >/dev/null 2>&1 || { echo >&2 "Error: minimap2 not found."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "Error: samtools not found."; exit 1; }

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
# Step 2: Alignment with minimap2
########################################
echo "Step 2: Aligning reads with minimap2..."

for i in {1..22}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    BAM_PREFIX="aligned_bams/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}"
    
    echo "Aligning ${FASTQ_FILTERED}..."
    minimap2 -ax map-ont \
             -t $THREADS \
             --MD \
             -Y \
             "$REFERENCE_GENOME" \
             "$FASTQ_FILTERED" | \
    samtools sort -@ $THREADS -o "${BAM_PREFIX}.aligned.bam"
    
    samtools index -@ $THREADS "${BAM_PREFIX}.aligned.bam"
    
    # Generate alignment stats
    samtools flagstat "${BAM_PREFIX}.aligned.bam" > "${BAM_PREFIX}.flagstat"
    samtools stats "${BAM_PREFIX}.aligned.bam" > "${BAM_PREFIX}.stats"
done

step_complete "2" "Success!!"

echo ""
echo "========================================"
echo "Output files: aligned_bams/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.aligned.bam"
echo "========================================"

########################################
# Step 3: Consensus generation with samtools
########################################
echo "Step 3: Generating consensus sequences..."

for i in {1..22}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BAM="aligned_bams/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.aligned.bam"
    CONSENSUS="consensus_sequences/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fasta"
    
    echo "Processing ${BAM}..."
    
    # Verify BAM is aligned to this reference
    if ! samtools view -H "$BAM" | grep -q "$(head -n1 "$REFERENCE_GENOME" | tr -d '>')"; then
        echo "ERROR: BAM not aligned to specified reference" >&2
        exit 1
    fi
    
    samtools consensus -a --show-ins no "$BAM" -o "$CONSENSUS"
    
    if [ ! -s "$CONSENSUS" ]; then
        echo "ERROR: Empty consensus for barcode ${BARCODE_PADDED}" >&2
        exit 1
    fi
done

step_complete "3" "Success!!"

echo ""
echo "========================================"
echo "Output files: consensus_sequences/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fasta"
echo "========================================"

########################################
# Step 4: Select RSVA OR RSVB consensus per sample
########################################
echo "Step 4: Select highest-quality RSVA/RSVB consensus per sample..."

> consensus_sequences/RSVA_consensus.fasta
> consensus_sequences/RSVB_consensus.fasta

for i in {1..22}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    INPUT_FASTA="consensus_sequences/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fasta"
    
    echo "Processing ${INPUT_FASTA}..."
    
    # Count non-N bases for each virus in the sample
    RSVA_COUNT=$(awk '/^>RSVA/{getline; seq=""; while ((getline line)>0) {if(line~/^>/){break}; seq=seq line}; gsub(/[^ACGT]/, "", seq); print length(seq)}' "$INPUT_FASTA")
    RSVB_COUNT=$(awk '/^>RSVB/{getline; seq=""; while ((getline line)>0) {if(line~/^>/){break}; seq=seq line}; gsub(/[^ACGT]/, "", seq); print length(seq)}' "$INPUT_FASTA")
    
    # Select the version with more ACGT bases (fewer Ns)
    if [ "$RSVA_COUNT" -ge "$RSVB_COUNT" ]; then
        # Use RSVA version (multi-line aware)
        awk -v bar="barcode${BARCODE_PADDED}" '
            /^>RSVA/ { 
                print $0 "|" bar
                while ((getline line) > 0) {
                    if (line ~ /^>/) exit
                    print line
                }
            }
        ' "$INPUT_FASTA" >> consensus_sequences/RSVA_consensus.fasta
    else
        # Use RSVB version (multi-line aware)
        awk -v bar="barcode${BARCODE_PADDED}" '
            /^>RSVB/ { 
                print $0 "|" bar
                while ((getline line) > 0) {
                    if (line ~ /^>/) exit
                    print line
                }
            }
        ' "$INPUT_FASTA" >> consensus_sequences/RSVB_consensus.fasta
    fi
done

RSVA_COUNT=$(grep -c "^>RSVA" consensus_sequences/RSVA_consensus.fasta)
RSVB_COUNT=$(grep -c "^>RSVB" consensus_sequences/RSVB_consensus.fasta)

step_complete "4" "Success!!: Selected RSVA ($RSVA_COUNT samples) and RSVB ($RSVB_COUNT samples)"

echo ""
echo "========================================"
echo "Output files:"
echo "- RSVA consensus: consensus_sequences/RSVA_consensus.fasta"
echo "- RSVB consensus: consensus_sequences/RSVB_consensus.fasta"
echo "========================================"



