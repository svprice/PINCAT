#!/bin/bash

# INPUTS:
# - raw partek file
# - output directory
# - classification


GROUP=$1
OUTPUT_DIR=$2
CLASS=$3

python -c 'import sys, preprocess_functions as pf; \
    pf.partekSpecificPreprocessing(sys.argv[1],sys.argv[2],sys.argv[3])' \
    "$GROUP" "$OUTPUT_DIR" "$CLASS"


python -c 'import sys, preprocess_functions as pf; \
    pf.downloadEnsemblMapData(sys.argv[1],sys.argv[2])' \
    "$GROUP" "$OUTPUT_DIR"
    
python -c 'import sys, preprocess_functions as pf; \
    pf.mergeToCreateEnsemblMap(sys.argv[1],sys.argv[2])' \
    "$GROUP" "$OUTPUT_DIR"

python -c 'import sys, preprocess_functions as pf; \
    pf.applyAllEnsemblMaps(sys.argv[1],sys.argv[2])' \
    "$GROUP" "$OUTPUT_DIR"

Rscript interactions_from_ensembl.R "$GROUP" "$OUTPUT_DIR"

python -c 'import sys, preprocess_functions as pf; \
    pf.getNetworkData(sys.argv[1],sys.argv[2])' \
    "$GROUP" "$OUTPUT_DIR"
