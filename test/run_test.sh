#!/bin/bash
set -e

# Quick test script for sc-sra-nf pipeline
# Usage: ./run_test.sh [star_index_path] OR [cellranger_ref_path]

echo "=== sc-sra-nf Pipeline Test ==="
echo

# Get the script's directory and move to project root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# Check if test file exists
if [ ! -f "test/SraRunTable.csv" ]; then
    echo "Error: test/SraRunTable.csv not found."
    exit 1
fi

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <star_index_path|cellranger_ref_path>"
    echo
    echo "Examples:"
    echo "  ./run_test.sh /path/to/star_index"
    echo "  ./run_test.sh /path/to/refdata-gex-GRCh38-2024-A"
    exit 1
fi

REF_PATH=$1
OUTDIR="test_results_$(date +%Y%m%d_%H%M%S)"

# Detect reference type
if [ -d "$REF_PATH/SA" ] || [ -f "$REF_PATH/SAindex" ]; then
    REF_TYPE="star"
    echo "Detected STAR index at: $REF_PATH"
    QUANT_PARAM="--star_index $REF_PATH"
elif [ -d "$REF_PATH/fasta" ] && [ -d "$REF_PATH/genes" ]; then
    REF_TYPE="cellranger"
    echo "Detected CellRanger reference at: $REF_PATH"
    QUANT_PARAM="--cellranger_ref $REF_PATH"
else
    echo "Error: Could not detect reference type."
    echo "Expected STAR index (with SA/SAindex files) or CellRanger reference (with fasta/ and genes/ dirs)"
    exit 1
fi

echo "Output directory: $OUTDIR"
echo

# Run pipeline
echo "Starting pipeline..."
nextflow run main.nf \
    --sra_table test/SraRunTable.csv \
    $QUANT_PARAM \
    --outdir "test/$OUTDIR" \
    -c test/test.config \
    -resume

echo
echo "=== Test Complete ==="
echo "Results saved to: $OUTDIR"
echo
echo "Check results:"
echo "  ls -lh $OUTDIR/"
echo "  cat $OUTDIR/summary/all_samples_summary.txt"
