use clap::Parser;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Arc;

/// Rust implementation of GoTChAo pipeline (Aggregated Output)
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(long)]
    barcode_fastq: String,
    #[arg(long)]
    sequence_fastq: String,
    #[arg(long)]
    whitelist: String,
    #[arg(long)]
    primer: String,
    #[arg(long)]
    ref_codon: String,
    #[arg(long)]
    mut_codon: String,
    #[arg(long)]
    mut_start: usize, // 1-based
    #[arg(long)]
    mut_end: usize,   // 1-based
    #[arg(long)]
    max_mismatches: u8,
    #[arg(long)]
    out: String,
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    rc_barcode: bool,
}

struct Counts {
    wt: u32,
    mut_: u32,
}

const COMPLEMENT: [u8; 256] = {
    let mut c = [0u8; 256];
    let mut i = 0;
    while i < 256 { c[i] = i as u8; i += 1; }
    c[b'A' as usize] = b'T'; c[b'a' as usize] = b'T';
    c[b'C' as usize] = b'G'; c[b'c' as usize] = b'G';
    c[b'G' as usize] = b'C'; c[b'g' as usize] = b'C';
    c[b'T' as usize] = b'A'; c[b't' as usize] = b'A';
    c[b'N' as usize] = b'N'; c[b'n' as usize] = b'N';
    c
};

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| COMPLEMENT[b as usize]).collect()
}

fn load_whitelist(path: &str) -> HashSet<Vec<u8>> {
    eprintln!("=== Loading Whitelist: {} ===", path);
    let file = File::open(path).expect("Cannot open whitelist");
    let reader = BufReader::new(file);
    let mut set = HashSet::new();
    for line in reader.lines() {
        let l = line.expect("Error reading line");
        if l.starts_with("Barcode") { continue; }
        let bc = l.split(',').next().unwrap().split('-').next().unwrap().trim();
        set.insert(bc.as_bytes().to_vec());
    }
    eprintln!("--> Loaded {} barcodes.", set.len());
    set
}

fn is_primed(seq: &[u8], primer: &[u8], max_miss: u8) -> bool {
    let mut miss = 0;
    if seq.len() < primer.len() { return false; }
    for i in 0..primer.len() {
        if seq[i] != primer[i] {
            miss += 1;
            if miss > max_miss { return false; }
        }
    }
    true
}

fn main() {
    let args = Args::parse();
    let whitelist = Arc::new(load_whitelist(&args.whitelist));
    
    let primer = args.primer.as_bytes().to_vec();
    let ref_codon = args.ref_codon.as_bytes().to_vec();
    let mut_codon = args.mut_codon.as_bytes().to_vec();
    let m_start = args.mut_start - 1;
    let m_end = args.mut_end; 

    eprintln!("=== Processing Reads (Rust Core: Aggregation Mode) ===");
    if args.rc_barcode { eprintln!("*** RC Barcodes: ON ***"); }

    let mut reader_r1 = parse_fastx_file(&args.barcode_fastq).expect("Invalid R1");
    let mut reader_r2 = parse_fastx_file(&args.sequence_fastq).expect("Invalid R2");

    let chunk_size = 10_000_000;
    let mut batch_data = Vec::with_capacity(chunk_size);
    let mut global_counts: HashMap<String, Counts> = HashMap::new();
    let mut total_scanned = 0;
    let mut total_matched = 0;

    loop {
        let r1_rec = reader_r1.next();
        let r2_rec = reader_r2.next();
        if r1_rec.is_none() || r2_rec.is_none() { break; }

        let r1 = r1_rec.unwrap().expect("R1 parse error");
        let r2 = r2_rec.unwrap().expect("R2 parse error");
        
        total_scanned += 1;

        // --- PRINTING LOGIC FIXED HERE ---
        // Only prints if exactly divisible by 20,000,000
        if total_scanned % 20_000_000 == 0 {
            eprintln!("   Scanned: {} | Matched: {}", total_scanned, total_matched);
        }

        let seq_r1 = r1.seq();
        if seq_r1.len() < 16 { continue; }
        let raw_bc = &seq_r1[0..16];
        let final_bc = if args.rc_barcode { reverse_complement(raw_bc) } else { raw_bc.to_vec() };

        if whitelist.contains(&final_bc) {
            total_matched += 1;
            batch_data.push((String::from_utf8_lossy(&final_bc).to_string(), r2.seq().to_vec()));
        }

        if batch_data.len() >= chunk_size {
            let results = process_batch(&batch_data, &primer, &ref_codon, &mut_codon, m_start, m_end, args.max_mismatches);
            for (bc, is_wt, is_mut) in results {
                let entry = global_counts.entry(bc).or_insert(Counts { wt: 0, mut_: 0 });
                entry.wt += is_wt as u32;
                entry.mut_ += is_mut as u32;
            }
            batch_data.clear();
        }
    }

    if !batch_data.is_empty() {
        let results = process_batch(&batch_data, &primer, &ref_codon, &mut_codon, m_start, m_end, args.max_mismatches);
        for (bc, is_wt, is_mut) in results {
            let entry = global_counts.entry(bc).or_insert(Counts { wt: 0, mut_: 0 });
            entry.wt += is_wt as u32;
            entry.mut_ += is_mut as u32;
        }
    }

    eprintln!("=== Writing Output to CSV ===");
    
    // --- NEW: Create Directory if it doesn't exist ---
    let path = std::path::Path::new(&args.out);
    if let Some(parent) = path.parent() {
        if !parent.exists() {
            std::fs::create_dir_all(parent).expect("Cannot create output directory");
        }
    }
    // -------------------------------------------------
    let mut writer = csv::Writer::from_path(&args.out).expect("Cannot create output file");
    writer.write_record(&["Barcode", "WT_Count", "MUT_Count", "log10_WT", "log10_MUT"]).unwrap();

    for (bc, counts) in global_counts {
        let log_wt = (counts.wt as f64 + 1.0).log10();
        let log_mut = (counts.mut_ as f64 + 1.0).log10();
        writer.write_record(&[bc, counts.wt.to_string(), counts.mut_.to_string(), format!("{:.4}", log_wt), format!("{:.4}", log_mut)]).unwrap();
    }
    eprintln!("=== DONE. Total Scanned: {} | Total Matched: {} ===", total_scanned, total_matched);
}

fn process_batch(
    batch: &[(String, Vec<u8>)],
    primer: &[u8], ref_c: &[u8], mut_c: &[u8],
    m_start: usize, m_end: usize, max_miss: u8
) -> Vec<(String, u8, u8)> {
    batch.par_iter().filter_map(|(bc, seq)| {
        if seq.len() < m_end { return None; }
        if !is_primed(seq, primer, max_miss) { return None; }
        let obs_codon = &seq[m_start..m_end];
        let is_wt = if obs_codon == ref_c { 1 } else { 0 };
        let is_mut = if obs_codon == mut_c { 1 } else { 0 };
        if is_wt == 0 && is_mut == 0 { return None; }
        Some((bc.clone(), is_wt, is_mut))
    }).collect()
}