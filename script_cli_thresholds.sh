#!/bin/bash
THREADS=2

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
# Step 5.5: Nextclade Analysis with all QC checks disabled
########################################
echo "Step 5.5: Running Nextclade with all QC checks disabled..."

# Make directory for Nextclade results
mkdir -p trees/nextclade_results

# Define the RSV types and their corresponding datasets
declare -A rsv_types=(
  ["RSVA"]="rsv_a"
  ["RSVB"]="rsv_b"
  ["RSVAD"]="rsv_a"
  ["RSVBD"]="rsv_b"
)

# Create a temporary pathogen.json file with default RSV QC parameters
QC_CONFIG=$(mktemp)
cat <<EOF > "$QC_CONFIG"
{
  "schemaVersion": "1.10.0",
  "metadata": {
    "title": "RSV Analysis",
    "description": "RSV analysis with default QC parameters"
  },
  "qc": {
    "privateMutations": {
      "enabled": true,
      "typical": 10,
      "cutoff": 20
    },
    "missingData": {
      "enabled": true,
      "missingDataThreshold": 300
    },
    "snpClusters": {
      "enabled": true,
      "windowSize": 50,
      "clusterCutOff": 5
    },
    "mixedSites": {
      "enabled": true,
      "mixedSitesThreshold": 5
    },
    "frameShifts": {
      "enabled": false,
      "frameShiftsThreshold": 0
    }
  },
  "geneMap": [],
  "primers": []
}
EOF

# Process each RSV type
for rsv_type in "${!rsv_types[@]}"; do
  input_file="trees/pooled_consensus/pooled_${rsv_type}.fasta"
  output_prefix="trees/nextclade_results/pooled_${rsv_type}"
  dataset="${rsv_types[$rsv_type]}"
  
  if [[ -f "$input_file" ]]; then
    echo "Processing $rsv_type with dataset $dataset"
    
    nextclade run \
      "$input_file" \
      --output-csv "${output_prefix}_nextclade.csv" \
      --output-json "${output_prefix}_nextclade.json" \
      --output-fasta "${output_prefix}_aligned.fasta" \
      --dataset-name "$dataset" \
      --input-pathogen-json "$QC_CONFIG" \
      --verbose
  else
    echo "Input file $input_file not found, skipping $rsv_type Nextclade analysis"
  fi
done

# Clean up temporary config file
rm "$QC_CONFIG"

step_complete "5.5" "Nextclade analysis complete with all QC checks disabled"

END_TIME=$(date +%s)
echo "Total pipeline runtime: $((END_TIME - START_TIME)) seconds"