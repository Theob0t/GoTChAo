# GoT-ChAo ⚡️🦀
### Genotyping of Targeted loci with single-cell Chromatin Accessibility (GoT-ChA) — Accelerated & Optimized

![CI Status](https://github.com/theob0t/GoTChAo/actions/workflows/ci.yml/badge.svg)
![Rust](https://img.shields.io/badge/Core-Rust-orange?style=flat-square&logo=rust)
![Python](https://img.shields.io/badge/Analysis-Python-blue?style=flat-square&logo=python)
![Docker](https://img.shields.io/badge/Container-Docker-2496ED?style=flat-square&logo=docker)
![Apptainer](https://img.shields.io/badge/HPC-Apptainer%20Compatible-blueviolet?style=flat-square)

```text
         ______     _______     ____          ___          
        / ____|  __|__   __|  / ____| |      /   \    __   
        | |  __  /  \  | |    | |     |__   /  ^  \  /  \  
        | | |_ || () | | |  - | |     |  \ /  /_\  \| () | 
        | |__| | \__/  | |    | |____ |  |/  _____  \\__/  
        \______|       |_|    \______ |  /__/     \__\     

              GoT...Ciao! (It's already done). 🤌
```

## 📖 Overview

**GoT-ChAo** is a high-performance reimplementation of the [Landau Lab GoTChA pipeline](https://github.com/landau-lab/Gotcha). It is designed to perform precise genotyping of single-cell RNA-seq data (10x Genomics) to distinguish between Wild Type (WT) and Mutant (MUT) cells based on specific primers and gene targets.

### 🚀 Why GoT-ChAo?
The original pipeline was written in **R/Python**. While accurate, it faced performance bottlenecks when processing massive datasets (100M+ reads). **GoT-ChAo** solves this by using a hybrid architecture that reduces runtime from **hours to minutes**:

1.  **Rust Core:** Handles heavy FASTQ parsing, barcode matching, and primer searching using **SIMD acceleration** and **Rayon parallelism**.
2.  **Optimized Python:** Re-engineers statistical genotyping with vectorized NumPy operations and analytical moments for noise correction.

| Feature | Original GoT-ChA | **GoT-ChAo** |
| :--- | :--- | :--- |
| **Language** | R / Python | **Rust** + Python |
| **Speed (60M Reads)** | Hours | **< 3 Minutes** |
| **Memory Usage** | High (Loads all reads) | **Low** (Streaming chunks) |
| **Deployment** | Conda Env | **Docker / Apptainer** |

---

## 🛠 Prerequisites

GoT-ChAo is containerized. You do not need Rust or Python installed. Choose your environment:

### 🅰️ For HPC Users (Clusters/Supercomputers)
If you are on a shared cluster (e.g., Slurm, LSF), you likely cannot run Docker. Use **Apptainer** (formerly Singularity).

*   **Requirements:** `apptainer` or `singularity` loaded.
*   **NYGC Users:** Load `squashfuse-default` via Conda or Modules before building.

### 🅱️ For Docker Users (Local Laptop/Cloud)
If you have root access or are running locally.

*   **Requirements:** [Docker Desktop](https://www.docker.com/products/docker-desktop/) or Docker Engine.

---

## 📦 Installation

We automatically build and verify images on GitHub Actions.

### Method 1: HPC / Apptainer (Recommended for Clusters)
You can pull the Docker image directly into a Singularity/Apptainer file (`.sif`).

```bash
# This converts the Docker image to a SIF file automatically
apptainer build gotchao.sif docker://ghcr.io/theob0t/gotchao:latest
```
*You now have `gotchao.sif` (approx 500MB) ready to run.*

### Method 2: Docker Pull
```bash
docker pull ghcr.io/theob0t/gotchao:latest
```

---

## 🏃‍♂️ Usage

### 1. Running on HPC (Apptainer)
The `.sif` file is a standalone executable. We use `-B` to bind your data folders so the container can see them.

```bash
# Example: Running on a cluster
# Note: We use 'apptainer run' to explicitly pass the -B bind flag
apptainer run -B /gpfs ./gotchao.sif \
         --barcode_fastq_path "/data/tests/data/test_barcode_R2.fastq.gz" \
         --sequence_fastq_path "/data/tests/data/test_seq_R1.fastq.gz" \
         --whitelist_path "/data/tests/data/test_singlecell.csv" \
         --primer_sequence CCTAGCCTGCCTCAGGAAACTGTGGATCAGGAACCCAAGGATCAGAA \
         --ref_codon GAG --mutation_codon GGA \
         --mutation_start 48 --mutation_end 50 --max_mismatches_primer 3 \
         --out "/data/tests/output_test_run"
```

### 2. Running with Docker
You must mount your current directory (volume) so results persist after the container stops.

```bash
# Mount current directory $(pwd) to /data inside container
docker run --rm -v $(pwd):/data ghcr.io/theob0t/gotchao:latest \
  --barcode_fastq_path /data/R2_001.fastq.gz \
  --sequence_fastq_path /data/R1_001.fastq.gz \
  --whitelist_path /data/singlecell.csv \
  --primer_sequence CCTAGCCTGCCTCAGGAAACTGTGGATCAGGAACCCAAGGATCAGAA \
  --ref_codon GAG \
  --mutation_codon GGA \
  --mutation_start 48 \
  --mutation_end 50 \
  --max_mismatches_primer 3 \
  --out /data/results_folder
```

---

## 🔧 Arguments Explained

| Argument | Description |
| :--- | :--- |
| `--barcode_fastq_path` | Path to the **Cell Barcode** FASTQ (usually R1 for 10x, but check your kit). |
| `--sequence_fastq_path` | Path to the **Biological Sequence** FASTQ (usually R2). |
| `--whitelist_path` | Path to CellRanger `singlecell.csv` (barcodes to analyze). |
| `--primer_sequence` | The anchor sequence to search for. |
| `--ref_codon` | Wild Type codon (e.g., `GAG`). |
| `--mutation_codon` | Mutant codon (e.g., `GGA`). |
| `--mutation_start` | **1-based** start index of the codon relative to the read. |
| `--mutation_end` | **1-based** end index of the codon. |
| `--out` | Directory for results. |
| `--no_rc` | Disable Reverse Complementing of barcodes (Default: ON). |

---

## 💻 Performance & Hardware

GoT-ChAo is extremely efficient compared to pure Python/R implementations.

*   **CPU:** Scales linearly. 4-8 cores is the sweet spot.
*   **RAM:** 
    *   **Rust Phase:** < 1GB (Streaming architecture).
    *   **Python Phase:** Depends on cell count (~4GB for 10k cells).
*   **Disk:** Read speed is the main bottleneck (~700k reads/sec on GPFS).

---

## 📊 Outputs

1.  `{SAMPLE}_counts.csv`: Raw aggregated counts (Rust output).
2.  `{SAMPLE}_genotype_labels.csv`: **Final Result.** Genotypes (`WT`, `MUT`, `HET`, `NA`) and confidence scores.
3.  `*_kde_mixture.pdf`: QC plot of noise distribution.
4.  `cluster_genotype.pdf`: Visualization of cell clusters.

---

## 📝 Attribution

*   **Original Concept:** [Landau Lab GoTChA](https://github.com/landau-lab/Gotcha)
```
