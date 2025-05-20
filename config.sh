#!/bin/bash

# Sample Prefix
SAMPLE_PREFIX=""

# Batch Name
BATCH=""

# Sample Name
declare -A SAMPLE_IDS=(
    ["barcode01"]="PMT"
    ["barcode02"]="PMT"
    ["barcode03"]="PMT"
    ["barcode04"]="PMT"
    ["barcode05"]="PMT"
    ["barcode06"]="PMT"
    ["barcode07"]="PMT"
    ["barcode08"]="PMT"
    ["barcode09"]="PMT"
    ["barcode10"]="PMT"
    ["barcode11"]="PMT"
    ["barcode12"]="PMT"
    ["barcode13"]="PMT"
    ["barcode14"]="PMT"
    ["barcode15"]="PMT"
    ["barcode16"]="PMT"
    ["barcode17"]="PMT"
    ["barcode18"]="PMT"
    ["barcode19"]="PMT"
    ["barcode20"]="PMT"
    ["barcode21"]="PMT"
    ["barcode22"]="PMT"
    ["barcode23"]="PMT"
    ["barcode24"]="PMT"
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
