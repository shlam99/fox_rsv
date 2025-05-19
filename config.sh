#!/bin/bash

# Sample Prefix
SAMPLE_PREFIX="b8140349-847c-4e49-a1cc-18f7e4f8d2be_SQK-NBD114-24"

# Batch Name
BATCH="20250425_1E"

# Sample Name
declare -A SAMPLE_IDS=(
    ["barcode01"]="PMT1111"
    ["barcode02"]="PMT2222"
    ["barcode03"]="PMT3333"
    ["barcode04"]="PMT3333"
    ["barcode05"]="PMT3333"
    ["barcode06"]="PMT3333"
    ["barcode07"]="PMT3333"
    ["barcode08"]="PMT3333"
    ["barcode09"]="PMT3333"
    ["barcode10"]="PMT3333"
    ["barcode11"]="PMT3333"
    ["barcode12"]="PMT3333"
    ["barcode14"]="PMT2222"
    ["barcode15"]="PMT2222"
    ["barcode16"]="PMT2222"
    ["barcode17"]="PMT2222"
    ["barcode18"]="PMT2222"
    ["barcode19"]="PMT2222"
    ["barcode20"]="PMT2222"
    ["barcode21"]="PMT2222"
    ["barcode22"]="PMT2222"
    ["barcode23"]="PMT2222"
    ["barcode24"]="PMT2222"
)

# Reference Genome
REFERENCE_GENOME=""

# Filtering Parameters
THREADS=4
MIN_QUALITY=7
MIN_LENGTH=500
MAX_LENGTH=50000
TARGET_BASES=100000000
KEEP_PERCENT=90
