#!/bin/bash

# INPUTS:
# - raw partek file
# - output directory
# - classification

RAW=$1
OUTPUT_DIR=$2
CLASS=$3

NEW_VAR=$(python -c 'import sys, preprocess_functions as pf; \
    pf.madeUpFunction(sys.argv[1],sys.argv[2],sys.argv[3])' \
    "$RAW" "$OUTPUT_DIR" "$CLASS")

echo $NEW_VAR
