
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
 
    # Add to appropriate consensus file with barcode prefix (Step 3 format) 
    echo "Adding ${BARCODE_ID} to ${RSV_TYPE} consensus..." 
    sed "s/^>/>${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}|/" "$CONSENSUS_SOURCE" >> "irma_consensus/${RSV_TYPE}_consensus.fasta" 
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
    SAMPLE_ID="${SAMPLE_IDS[$BARCODE_ID]}" 
    if [ -z "$SAMPLE_ID" ]; then 
        SAMPLE_ID="${BARCODE_ID}" 
    fi 
     
    # Add to batch file with sample ID (Step 4 format) 
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
