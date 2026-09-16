use clap::Parser;
use cpg_likelihoods::{append_pl_columns, estimate_error_c2n_from_bedmethyl};

#[derive(Parser, Debug)]
#[command(
    author = "Søren Blikdal Hanssen",
    version = "0.1.0",
    about = "Compute e_c2n and append PL_hom_cpg/PL_het_cpg/PL_non_cpg as cols 19-21"
)]
struct Args {
    /// Input bedMethyl file
    path: String,

    /// Output file with PL columns appended
    #[arg(long)]
    out: String,

    /// Minimum coverage for e_c2n estimation
    #[arg(long, default_value_t = 10)]
    min_cov: u64,

    /// Maximum coverage for e_c2n estimation
    #[arg(long, default_value_t = 100)]
    max_cov: u64,

    /// Minimum ratio cov/(cov+d) for e_c2n estimation
    #[arg(long, default_value_t = 0.8)]
    min_ratio: f64,

    /// n parameter in likelihood model
    #[arg(long, default_value_t = 2.0)]
    n: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let e_c2n = estimate_error_c2n_from_bedmethyl(
        &args.path,
        args.min_cov,
        args.max_cov,
        args.min_ratio,
    )?;

    eprintln!("e_c2n = {:.10}", e_c2n);

    append_pl_columns(&args.path, &args.out, e_c2n, args.n)?;

    eprintln!("Wrote: {}", args.out);
    Ok(())
}