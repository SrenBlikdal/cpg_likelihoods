use clap::Parser;
use cpg_likelihoods::{classify_from_pl, open_bedmethyl_reader, LocusClass};
use std::io::BufRead;

#[derive(Parser, Debug)]
#[command(
    author = "Søren Blikdal Hanssen",
    version = "0.1.0",
    about = "Locus classification summary from PL cols 19-21"
)]
struct Args {
    /// Input bedMethyl file with PL columns in 19-21 (.bed, .bed.gz, or .bed.bgz)
    path: String,

    /// Minimum GQ threshold
    #[arg(long, default_value_t = 20)]
    min_gq: u32,
}

#[derive(Default)]
struct Counts {
    total: u64,
    hom: u64,
    het: u64,
    non: u64,
    unclassifiedcpg: u64,
    unclassified: u64,
}

fn pct(n: u64, d: u64) -> f64 {
    if d == 0 { 0.0 } else { (n as f64) * 100.0 / (d as f64) }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1-based cols 19,20,21 => 0-based 18,19,20
    const PL_HOM_IDX: usize = 18;
    const PL_HET_IDX: usize = 19;
    const PL_NON_IDX: usize = 20;
    const MIN_FIELDS: usize = 21;

    let args = Args::parse();

    let reader = open_bedmethyl_reader(&args.path)?;

    let mut c = Counts::default();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < MIN_FIELDS {
            continue;
        }

        let pl_hom: u32 = match fields[PL_HOM_IDX].parse() {
            Ok(v) => v,
            Err(_) => continue, // header-like row
        };
        let pl_het: u32 = match fields[PL_HET_IDX].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let pl_non: u32 = match fields[PL_NON_IDX].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        c.total += 1;
        match classify_from_pl(pl_hom, pl_het, pl_non, args.min_gq) {
            LocusClass::HomCpg => c.hom += 1,
            LocusClass::HetCpg => c.het += 1,
            LocusClass::NonCpg => c.non += 1,
            LocusClass::UnclassifiedCpg => c.unclassifiedcpg += 1,
            LocusClass::Unclassified => c.unclassified += 1,
        }
    }

    let hom_het = c.hom + c.het + c.unclassifiedcpg;
    let all_unclassified = c.unclassified + c.unclassifiedcpg;

    println!("Input: {}", args.path);
    println!("min_GQ: {}", args.min_gq);
    println!("Total loci: {}", c.total);
    println!();
    println!("Classification with separate homozygous and heterozygous CpG states:");  
    println!("hom_cpg         {:>12}  {:>7.3}%", c.hom, pct(c.hom, c.total));
    println!("het_cpg         {:>12}  {:>7.3}%", c.het, pct(c.het, c.total));
    println!("non_cpg         {:>12}  {:>7.3}%", c.non, pct(c.non, c.total));
    println!("unclassified    {:>12}  {:>7.3}%", all_unclassified, pct(all_unclassified, c.total));
    println!();
    println!("Classification with collapsed homozygous and heterozygous CpG states:"); 
    println!("hom/het_cpg     {:>12}  {:>7.3}%", hom_het, pct(hom_het, c.total));
    println!("non_cpg         {:>12}  {:>7.3}%", c.non, pct(c.non, c.total));
    println!("unclassified    {:>12}  {:>7.3}%", c.unclassified, pct(c.unclassified, c.total));

    Ok(())
}