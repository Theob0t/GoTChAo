#!/bin/bash

# 1. robustly find the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# 2. Define paths relative to the repository
# This allows it to work on any computer (yours, colleague's, GitHub)
SIF_IMAGE="${REPO_ROOT}/gotchao.sif"
DATA_DIR="${SCRIPT_DIR}/data"
OUTPUT_DIR="${SCRIPT_DIR}/output_test_run"

# 3. Detect Container Runtime (Apptainer vs Singularity)
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
else
    CONTAINER_CMD="singularity"
fi

# 4. Handle Bind Mounts Smartly
# We always bind the current repo so it can write outputs
BIND_ARGS="-B ${REPO_ROOT}"

# If /gpfs exists (like on your HPC), bind it too. 
# If not (like on GitHub), skip it so the pipeline doesn't crash.
if [ -d "/gpfs" ]; then
    BIND_ARGS="${BIND_ARGS} -B /gpfs"
fi

echo "Starting Test Run..."
echo "Container: $CONTAINER_CMD"
echo "Image: $SIF_IMAGE"

mkdir -p "$OUTPUT_DIR"

# 5. THE RUN COMMAND
$CONTAINER_CMD run $BIND_ARGS "$SIF_IMAGE" \
  --barcode_fastq_path "${DATA_DIR}/test_barcode_R2.fastq.gz" \
  --sequence_fastq_path "${DATA_DIR}/test_seq_R1.fastq.gz" \
  --whitelist_path "${DATA_DIR}/test_whitelist.csv" \
  --primer_sequence CCTAGCCTGCCTCAGGAAACTGTGGATCAGGAACCCAAGGATCAGAA \
  --ref_codon GAG \
  --mutation_codon GGA \
  --mutation_start 48 \
  --mutation_end 50 \
  --max_mismatches_primer 3 \
  --out "$OUTPUT_DIR"

echo ""
echo "✅ Test run finished."
echo "   Check results in: $OUTPUT_DIR"