#!/bin/bash
set -e # Stop script if any command fails

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${SCRIPT_DIR}/data"
OUTPUT_DIR="${SCRIPT_DIR}/output_test_run"

# Clean previous run
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

echo "--- 1. Running Pipeline ---"

# Logic: If we are in CI (GitHub), use Docker. If on HPC, use Apptainer/Singularity.
if [ "$CI" = "true" ]; then
    # Running in GitHub Actions
    echo "Running in Docker Mode..."
    docker run --rm \
        -v "${REPO_ROOT}:/data" \
        gotchao-image \
        --barcode_fastq_path "/data/tests/data/R2_001.fastq.gz" \
        --sequence_fastq_path "/data/tests/data/R1_001.fastq.gz" \
        --whitelist_path "/data/tests/data/singlecell.csv" \
        --primer_sequence CCTAGCCTGCCTCAGGAAACTGTGGATCAGGAACCCAAGGATCAGAA \
        --ref_codon GAG --mutation_codon GGA \
        --mutation_start 48 --mutation_end 50 --max_mismatches_primer 3 \
        --out "/data/tests/output_test_run"
else
    # Running locally or HPC 
    echo "Please run via Apptainer/Singularity manually"
fi

echo "--- 2. Verifying Results ---"
# Compare the output we just made against the 'Gold Standard' committed in git
python3 tests/verify_results.py \
    "${OUTPUT_DIR}/YOUR_OUTPUT_FILE.csv" \
    "${DATA_DIR}/expected_output.csv"