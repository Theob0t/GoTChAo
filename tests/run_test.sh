#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status

# 1. Setup Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${SCRIPT_DIR}/data"
OUTPUT_DIR="${SCRIPT_DIR}/output_test_run"

# Gold Standard Files
EXPECTED_FILE="${DATA_DIR}/expected_output.csv"
RESULT_FILENAME="test_seq_R1_genotype_labels.csv"

# Clean previous run
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

echo "--- 1. Running Pipeline ---"

if [ "$CI" = "true" ]; then
    # ---------------------------------------------------------
    # GITHUB ACTIONS (DOCKER)
    # ---------------------------------------------------------
    echo "Environment: CI detected. Running in Docker Mode..."
    
    docker run --rm \
        -v "${REPO_ROOT}:/data" \
        gotchao-image \
        --barcode_fastq_path "/data/tests/data/test_barcode_R2.fastq.gz" \
        --sequence_fastq_path "/data/tests/data/test_seq_R1.fastq.gz" \
        --whitelist_path "/data/tests/data/test_singlecell.csv" \
        --primer_sequence CCTAGCCTGCCTCAGGAAACTGTGGATCAGGAACCCAAGGATCAGAA \
        --ref_codon GAG --mutation_codon GGA \
        --mutation_start 48 --mutation_end 50 --max_mismatches_primer 3 \
        --out "/data/tests/output_test_run"
    
    echo "--- 2. Verifying Results (Inside Docker) ---"
    if [ -f "$EXPECTED_FILE" ]; then
        # We override the ENTRYPOINT to run python directly on the verification script
        docker run --rm --entrypoint python3 \
            -v "${REPO_ROOT}:/data" \
            gotchao-image \
            "/data/tests/verify_results.py" \
            "/data/tests/output_test_run/${RESULT_FILENAME}" \
            "/data/tests/data/expected_output.csv"
    else
        echo "WARNING: No expected_output.csv found. Skipping."
    fi

else
    # ---------------------------------------------------------
    # HPC / LOCAL (APPTAINER)
    # ---------------------------------------------------------
    echo "Environment: Local/HPC detected."
    
    # 1. Locate the Image
    SIF_IMAGE="${REPO_ROOT}/gotchao.sif"
    
    if [ ! -f "$SIF_IMAGE" ]; then
        echo "ERROR: ${SIF_IMAGE} not found!"
        echo "Please build it first: apptainer build gotchao.sif docker://ghcr.io/theob0t/gotchao:latest"
        exit 1
    fi

    # 2. Detect Command
    if command -v apptainer &> /dev/null; then
        CONTAINER_CMD="apptainer"
    elif command -v singularity &> /dev/null; then
        CONTAINER_CMD="singularity"
    else
        echo "ERROR: Neither apptainer nor singularity found."
        exit 1
    fi

    echo "Using: $CONTAINER_CMD"
    echo "Image: $SIF_IMAGE"

    # 3. Run
    $CONTAINER_CMD run -B "${REPO_ROOT}" "$SIF_IMAGE" \
        --barcode_fastq_path "${DATA_DIR}/test_barcode_R2.fastq.gz" \
        --sequence_fastq_path "${DATA_DIR}/test_seq_R1.fastq.gz" \
        --whitelist_path "${DATA_DIR}/test_singlecell.csv" \
        --primer_sequence CCTAGCCTGCCTCAGGAAACTGTGGATCAGGAACCCAAGGATCAGAA \
        --ref_codon GAG --mutation_codon GGA \
        --mutation_start 48 --mutation_end 50 --max_mismatches_primer 3 \
        --out "${OUTPUT_DIR}"

    echo "--- 2. Verifying Results (Inside Container) ---"
    if [ -f "$EXPECTED_FILE" ]; then
        # Use 'exec' to run python3 inside the container
        $CONTAINER_CMD exec -B "${REPO_ROOT}" "$SIF_IMAGE" \
            python3 "${SCRIPT_DIR}/verify_results.py" \
            "${OUTPUT_DIR}/${RESULT_FILENAME}" \
            "$EXPECTED_FILE"
    else
        echo "WARNING: No "$EXPECTED_FILE" found. Skipping verification."
    fi
fi