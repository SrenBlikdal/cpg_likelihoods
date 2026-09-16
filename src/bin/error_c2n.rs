use clap::Parser;
use cpg_likelihoods::estimate_error_c2n_from_bedmethyl;

#[derive(Parser, Debug)]
#[command(
    author = "Søren Blikdal Hanssen",
    version = "0.1.0",
    about = "Estimate CpG->non-CpG error rate (e_c2n) from bedMethyl"
)]
struct Args {
    /// Input bedMethyl file (.bed, .bed.gz, or .bed.bgz)
    path: String,

    /// Minimum coverage (col5)
    #[arg(long, default_value_t = 10)]
    min_cov: u64,

    /// Maximum coverage (col5)
    #[arg(long, default_value_t = 100)]
    max_cov: u64,

    /// Minimum ratio cov/(cov+d)
    #[arg(long, default_value_t = 0.8)]
    min_ratio: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let e = estimate_error_c2n_from_bedmethyl(
        &args.path,
        args.min_cov,
        args.max_cov,
        args.min_ratio,
    )?;

    println!("{:.5}", e);
    Ok(())
}