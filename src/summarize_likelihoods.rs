use cpg_likelihoods::{BedMethylRecord, LocusClass, classify_from_pl, open_bedmethyl_reader};
use std::io::BufRead;
use std::path::Path;

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
    if d == 0 {
        0.0
    } else {
        (n as f64) * 100.0 / (d as f64)
    }
}

pub fn run(input: &Path, min_gq: u32) -> Result<(), Box<dyn std::error::Error>> {
    let reader = open_bedmethyl_reader(input)?;

    let mut c = Counts::default();

    for (index, line) in reader.lines().enumerate() {
        let context = |error| format!("{}:{}: {error}", input.display(), index + 1);
        let line = line.map_err(|error| context(error.to_string()))?;
        let Some(record) = BedMethylRecord::parse(&line).map_err(context)? else {
            continue;
        };
        let [pl_hom, pl_het, pl_non] = record.pl().map_err(context)?;

        c.total += 1;
        match classify_from_pl(pl_hom, pl_het, pl_non, min_gq) {
            LocusClass::HomCpg => c.hom += 1,
            LocusClass::HetCpg => c.het += 1,
            LocusClass::NonCpg => c.non += 1,
            LocusClass::UnclassifiedCpg => c.unclassifiedcpg += 1,
            LocusClass::Unclassified => c.unclassified += 1,
        }
    }

    let hom_het = c.hom + c.het + c.unclassifiedcpg;
    let all_unclassified = c.unclassified + c.unclassifiedcpg;

    println!("Input: {}", input.display());
    println!("min_GQ: {}", min_gq);
    println!("Total loci: {}", c.total);
    println!();
    println!("Classification with separate homozygous and heterozygous CpG states:");
    println!(
        "hom_cpg         {:>12}  {:>7.3}%",
        c.hom,
        pct(c.hom, c.total)
    );
    println!(
        "het_cpg         {:>12}  {:>7.3}%",
        c.het,
        pct(c.het, c.total)
    );
    println!(
        "non_cpg         {:>12}  {:>7.3}%",
        c.non,
        pct(c.non, c.total)
    );
    println!(
        "unclassified    {:>12}  {:>7.3}%",
        all_unclassified,
        pct(all_unclassified, c.total)
    );
    println!();
    println!("Classification with collapsed homozygous and heterozygous CpG states:");
    println!(
        "hom/het_cpg     {:>12}  {:>7.3}%",
        hom_het,
        pct(hom_het, c.total)
    );
    println!(
        "non_cpg         {:>12}  {:>7.3}%",
        c.non,
        pct(c.non, c.total)
    );
    println!(
        "unclassified    {:>12}  {:>7.3}%",
        c.unclassified,
        pct(c.unclassified, c.total)
    );

    Ok(())
}
