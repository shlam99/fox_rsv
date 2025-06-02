#!/bin/bash

# Sample Prefix
SAMPLE_PREFIX=""

# Batch Name
BATCH=""

# Sample Name
get_sample_id() {
    local barcode_id="$1"
    case "$barcode_id" in
        "barcode01") echo "PMT" ;;
        "barcode02") echo "PMT" ;;
        "barcode03") echo "PMT" ;;
        "barcode04") echo "PMT" ;;
        "barcode05") echo "PMT" ;;
        "barcode06") echo "PMT" ;;
        "barcode07") echo "PMT" ;;
        "barcode08") echo "PMT" ;;
        "barcode09") echo "PMT" ;;
        "barcode10") echo "PMT" ;;
        "barcode11") echo "PMT" ;;
        "barcode12") echo "PMT" ;;
        "barcode13") echo "PMT" ;;
        "barcode14") echo "PMT" ;;
        "barcode15") echo "PMT" ;;
        "barcode16") echo "PMT" ;;
        "barcode17") echo "PMT" ;;
        "barcode18") echo "PMT" ;;
        "barcode19") echo "PMT" ;;
        "barcode20") echo "PMT" ;;
        "barcode21") echo "PMT" ;;
        "barcode22") echo "PMT" ;;
        "barcode23") echo "PMT" ;;
        "barcode24") echo "PMT" ;;
        *) echo "$barcode_id" ;;  # Fallback to barcode ID if no match
    esac
}

# Reference Genome
REFERENCE_GENOME=""

# Filtering Parameters
THREADS=$(nproc --all)
MIN_QUALITY=7
MIN_LENGTH=500
MAX_LENGTH=50000
TARGET_BASES=100000000
KEEP_PERCENT=90
