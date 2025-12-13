# GoT-ChAo ⚡️🦀
### Genotyping of Targeted loci with single-cell Chromatin Accessibility (GoT-ChA) — Accelerated & Optimized

![Rust](https://img.shields.io/badge/Core-Rust-orange?style=flat-square&logo=rust)
![Python](https://img.shields.io/badge/Analysis-Python-blue?style=flat-square&logo=python)
![Container](https://img.shields.io/badge/Container-Apptainer%2FSingularity-blueviolet?style=flat-square)
![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)

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
The original pipeline was written in **R/Python**. While accurate, it faced performance bottlenecks when processing massive datasets (100M+ reads), taking hours to complete. **GoT-ChAo** solves this by using a hybrid architecture that reduces runtime to minutes:

1.  **Rust Core (`gotchao_core`):** Handles the heavy lifting (FASTQ parsing, barcode matching, primer search, and counting). It uses **SIMD acceleration**, **Zero-Copy memory management**, and **Rayon parallelism**.
2.  **Optimized Python Analytics:** Re-engineers the original statistical genotyping logic for speed:
    *   **Vectorization:** Replaces row-by-row iteration with instantaneous NumPy matrix operations.
    *   **Analytical Moments:** Calculates noise distribution statistics analytically rather than using slow numerical integration.
    *   **One-Pass KNN:** Optimizes the label propagation step to classify uncertain cells in a single pass using confident "anchor" cells.

| Feature | Original GoTChA | **GoT-ChAo** |
| :--- | :--- | :--- |
| **Language** | R / Python | **Rust** + Python |
| **Speed (60M Reads)** | Hours | **< 3 Minutes** |
| **Memory Usage** | High (Loads all reads) | **Low** (Streaming chunks) |
| **Parallelism** | Multiprocessing (Pickle overhead) | **True Multithreading** |
| **Deployment** | Conda Env | **Singularity / Apptainer** |

---

## ⚙️ Architecture

1.  **Input:** Raw FASTQ files (R1/R2) + CellRanger Whitelist (`singlecell.csv`).
2.  **Phase 1 (Rust):** 
    *   Streams FASTQ files in 1M read chunks.
    *   Matches Cell Barcodes (with optional Reverse Complement).
    *   Matches Primer sequences (allowing **3 mismatches** by default).
    *   Extracts the specific codon at the mutation site.
    *   Aggregates counts (WT vs MUT) per barcode in memory.
3.  **Phase 2 (Python):**
    *   Loads aggregated counts.
    *   Performs Noise Correction using Kernel Density Estimation (KDE).
    *   Classifies cells using optimized Quadrant gating and KNN refinement.
    *   Outputs final Genotype labels (`WT`, `MUT`, `HET`, `NA`).

---

## 🛠 Prerequisites

The only requirement to run this pipeline is **Apptainer** (formerly Singularity).

### For NYGC Users (New York Genome Center)
The NYGC HPC clusters require `squashfuse` to build/run containers.
```bash
# Activate the required module before running
conda activate squashfuse-default
```

### For General HPC Users
Ensure `apptainer` or `singularity` is installed and available in your path.
```bash
module load apptainer
# OR
module load singularity
```

---

## 📦 Installation

### Option 1: Quick Start (Recommended)
You do not need to install Rust or Python manually. Simply download the pre-built container from the releases page.

1.  Go to the [**Releases Page**](https://github.com/theob0t/GoTChAo/releases).
2.  Download the latest `gotchao.sif` file.
3.  You are ready to run.

### Option 2: Build from Source
If you want to modify the code or build it yourself:

```bash
git clone https://github.com/theob0t/GoTChAo.git
cd GoTChAo

# Build the image (takes ~2-3 minutes)
# --fakeroot allows building without sudo privileges
singularity build --fakeroot gotchao.sif Singularity.def
```

---

## 🏃‍♂️ Usage

To run GoT-ChAo, simply execute the container using `singularity run`.

**Note:** You must replace the `<PLACEHOLDERS>` below with the specific sequences and coordinates for your target gene.

```bash
# Note: -B /gpfs binds your filesystem so the container can see your data
singularity run -B /gpfs ./gotchao.sif \
  --barcode_fastq_path /path/to/R2_001.fastq.gz \
  --sequence_fastq_path /path/to/R1_001.fastq.gz \
  --whitelist_path /path/to/outs/singlecell.csv \
  --primer_sequence <YOUR_PRIMER_SEQUENCE_HERE> \
  --ref_codon <WT_CODON> \
  --mutation_codon <MUT_CODON> \
  --mutation_start <MUT_START_INDEX> \
  --mutation_end <MUT_END_INDEX> \
  --max_mismatches_primer 3 \
  --out /path/to/output_directory
```

### Arguments Explained
| Argument | Description |
| :--- | :--- |
| `--barcode_fastq_path` | Path to the file containing the **Cell Barcode** (usually R1 for standard 10x, but R2 for some custom libraries). |
| `--sequence_fastq_path` | Path to the file containing the **Biological Sequence** (usually R2 for standard 10x, but R1 for some). |
| `--whitelist_path` | Path to the CellRanger `singlecell.csv` output. |
| `--primer_sequence` | The anchor sequence to search for (allowing 3 mismatches by default). |
| `--ref_codon` | The Wild Type codon sequence (e.g., `GAG`). |
| `--mutation_codon` | The Mutant codon sequence (e.g., `GGA`). |
| `--mutation_start` | **1-based** start index of the codon relative to the read start. |
| `--mutation_end` | **1-based** end index of the codon. |
| `--out` | **Directory** where results will be saved. |
| `--no_rc` | (Optional) Use this flag to **disable** Reverse Complementing barcodes. *Default is ON.* |

---

## 🧪 Testing

We provide a micro-dataset (cherry-picked real reads) to verify the pipeline works.

1.  **Generate Test Data:** (If not already present)
    ```bash
    # Only needed if you have access to the source data
    python tests/make_test_data.py
    ```
2.  **Run Test:**
    ```bash
    ./tests/run_test.sh
    ```
    *   **Expected Time:** < 5 seconds.
    *   **Expected Output:** A CSV in `tests/data/results/` with properly labeled WT and MUT cells.

---

## 💻 Hardware Requirements

GoT-ChAo is extremely efficient.

*   **CPU:** 4+ Cores recommended (It scales linearly with cores).
*   **RAM:** 
    *   **Rust Phase:** Minimal (~500MB). It streams data.
    *   **Python Phase:** Depends on the number of cells. For 10k cells, < 4GB RAM. For 1M cells, ~16GB RAM is recommended.
*   **Disk:** Read speed is the main bottleneck. It can process ~700,000 reads/second on standard GPFS storage.

---

## 📊 Outputs

The output directory will contain:
1.  `{SAMPLE}_counts.csv`: Raw aggregated counts from the Rust core.
2.  `{SAMPLE}_genotype_labels.csv`: **Final Result.** Contains Genotype (`WT`, `MUT`, `HET`, `NA`), Confidence scores, and transformed counts.
3.  `*_kde_mixture.pdf`: QC plot showing the noise threshold determination.
4.  `cluster_genotype.pdf`: Scatter plot of the final clusters.

---

## 📝 Attribution

*   **Original Concept:** [Landau Lab GoTChA](https://github.com/landau-lab/Gotcha)
