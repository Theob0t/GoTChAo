import argparse
import subprocess
import sys
import os
import time



print(r'''
 ______     _______     ____          ___          
/ ____|  __|__   __|  / ____| |      /   \    __   
| |  __  /  \  | |    | |     |__   /  ^  \  /  \  
| | |_ || () | | |  - | |     |  \ /  /_\  \| () | 
| |__| | \__/  | |    | |____ |  |/  _____  \\__/  
\______|       |_|    \______ |  /__/     \__\     
''')


script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

try:
    from gotchao_labeling import GotchaLabeling
except ImportError:
    print("Warning: Could not import 'gotchao_labeling'. Genotyping step will fail.")

def main():
    pipeline_start = time.time()

    parser = argparse.ArgumentParser(description="GoTChAo Pipeline (Rust Core + Python Genotyping)")
    parser.add_argument('--barcode_fastq_path', required=True)
    parser.add_argument('--sequence_fastq_path', required=True)
    parser.add_argument('--whitelist_path', required=True)
    parser.add_argument('--primer_sequence', required=True)
    parser.add_argument('--ref_codon', required=True)
    parser.add_argument('--mutation_codon', required=True)
    parser.add_argument('--mutation_start', required=True)
    parser.add_argument('--mutation_end', required=True)
    parser.add_argument('--max_mismatches_primer', required=True)
    
    # UPDATED HELP TEXT
    parser.add_argument('--out', required=True, help="Output DIRECTORY (not file)")
    
    parser.add_argument('--no_rc', dest='rc_barcode', action='store_false')
    parser.set_defaults(rc_barcode=True)

    args = parser.parse_args()
    
    # --- SETUP PATHS ---
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    
    # Infer Sample ID from the Sequence FASTQ filename
    # Example: "CH40_LM250318...fastq.gz" -> "CH40_LM250318"
    fastq_name = os.path.basename(args.sequence_fastq_path)
    sample_id = fastq_name.split('.')[0]
    
    # Simplify ID if it's too long (split by first few underscores?)
    # For now, let's just use the prefix before first period
    
    counts_filename = f"{sample_id}_counts.csv"
    counts_path = os.path.join(out_dir, counts_filename)
    
    rust_bin = "/usr/local/bin/gotchao_core"
    
    print("\n=== STEP 1: Launching Rust Core (Counting) ===")
    print(f"   Output Counts File: {counts_path}")
    
    cmd = [
        rust_bin,
        "--barcode-fastq", args.barcode_fastq_path,
        "--sequence-fastq", args.sequence_fastq_path,
        "--whitelist", args.whitelist_path,
        "--primer", args.primer_sequence,
        "--ref-codon", args.ref_codon,
        "--mut-codon", args.mutation_codon,
        "--mut-start", str(args.mutation_start),
        "--mut-end", str(args.mutation_end),
        "--max-mismatches", str(args.max_mismatches_primer),
        "--out", counts_path, # Pass the FILE path to Rust
        f"--rc-barcode={'true' if args.rc_barcode else 'false'}"
    ]
    
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        print(f"Error: Rust binary not found at {rust_bin}. Running inside container?")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Rust binary failed with code {e.returncode}")
        sys.exit(e.returncode)

    # --- STEP 2: PYTHON GENOTYPING (Labeling) ---
    print("\n=== STEP 2: Launching Python Genotyping (Labeling) ===")
    print(f"   Working Directory: {out_dir}")
    print(f"   Input File: {counts_filename}")
    
    try:
        GotchaLabeling(
            path=out_dir,
            infile=counts_filename,
            sample_id=sample_id
        )
    except Exception as e:
        print(f"Genotyping Error: {e}")
        sys.exit(1)

    pipeline_end = time.time()
    elapsed = pipeline_end - pipeline_start
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)

    print("\n" + "="*40)
    print(f"✅ PIPELINE COMPLETE in {minutes}m {seconds}s")
    print(f"   Results are in: {out_dir}")
    print("="*40 + "\n")

if __name__ == "__main__":
    main()