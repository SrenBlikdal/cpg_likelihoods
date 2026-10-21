use clap::{Args, Parser, Subcommand};
use cpg_likelihoods::{
    FilterState, append_pl_columns, estimate_error_c2n_from_bedmethyl, filter_bedmethyl_by_state,
};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "cpg-likelihoods",
    version,
    about = "Infer CpG states from bedMethyl counts"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Estimate the CpG-to-non-CpG error rate
    #[command(name = "error_c2n")]
    ErrorC2n(ErrorC2nArgs),
    /// Append Phred-scaled likelihoods as columns 19-21
    Likelihoods(LikelihoodsArgs),
    /// Summarize CpG states from a PL-annotated bedMethyl file
    #[command(name = "summarize_likelihoods")]
    SummarizeLikelihoods(SummarizeLikelihoodsArgs),
    /// Filter a PL-annotated bedMethyl file by CpG state
    #[command(name = "filter_likelihoods")]
    FilterLikelihoods(FilterLikelihoodsArgs),
}

impl Command {
    pub fn run(self) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            Self::ErrorC2n(args) => args.run(),
            Self::Likelihoods(args) => args.run(),
            Self::SummarizeLikelihoods(args) => args.run(),
            Self::FilterLikelihoods(args) => args.run(),
        }
    }
}

#[derive(Debug, Args)]
pub struct ErrorC2nArgs {
    /// Input bedMethyl file (plain, gzip, or BGZF)
    pub input: PathBuf,
    #[command(flatten)]
    pub estimation: EstimationArgs,
}

impl ErrorC2nArgs {
    pub fn run(self) -> Result<(), Box<dyn std::error::Error>> {
        let e = estimate_error_c2n_from_bedmethyl(
            &self.input,
            self.estimation.min_cov,
            self.estimation.max_cov,
            self.estimation.min_ratio,
        )?;
        println!("{e:.10}");
        Ok(())
    }
}

#[derive(Debug, Args)]
pub struct LikelihoodsArgs {
    /// Input bedMethyl file (plain, gzip, or BGZF)
    pub input: PathBuf,
    /// Output file (.bed, .bed.gz, or .bed.bgz)
    #[arg(long)]
    pub out: PathBuf,
    /// Supply an error rate instead of estimating it
    #[arg(long)]
    pub e_c2n: Option<f64>,
    #[command(flatten)]
    pub estimation: EstimationArgs,
    /// Common likelihood scale factor (does not affect PLs)
    #[arg(long, default_value_t = 2.0)]
    pub n: f64,
}

impl LikelihoodsArgs {
    pub fn run(self) -> Result<(), Box<dyn std::error::Error>> {
        let e = match self.e_c2n {
            Some(e) => {
                if !(0.0..=1.0).contains(&e) {
                    return Err("e_c2n must be between 0 and 1".into());
                }
                eprintln!("Using manual e_c2n = {e:.10}");
                e
            }
            None => {
                let e = estimate_error_c2n_from_bedmethyl(
                    &self.input,
                    self.estimation.min_cov,
                    self.estimation.max_cov,
                    self.estimation.min_ratio,
                )?;
                eprintln!("Estimated e_c2n = {e:.10}");
                e
            }
        };
        append_pl_columns(&self.input, &self.out, e, self.n)?;
        eprintln!("Wrote: {}", self.out.display());
        Ok(())
    }
}

#[derive(Debug, Args)]
pub struct SummarizeLikelihoodsArgs {
    pub input: PathBuf,
    #[arg(long, default_value_t = 20)]
    pub min_gq: u32,
}

impl SummarizeLikelihoodsArgs {
    pub fn run(self) -> Result<(), Box<dyn std::error::Error>> {
        crate::summarize_likelihoods::run(&self.input, self.min_gq)
    }
}

#[derive(Debug, Args)]
pub struct FilterLikelihoodsArgs {
    #[arg(long)]
    pub input_bed: PathBuf,
    #[arg(long)]
    pub output_bed: PathBuf,
    /// CpG_hom, CpG_het, non-CpG, or CpG_hom/het
    #[arg(long)]
    pub state: FilterState,
    #[arg(long, default_value_t = 20)]
    pub min_gq: u32,
    /// Retain the three PL columns in the output
    #[arg(long)]
    pub keep_pl: bool,
}

impl FilterLikelihoodsArgs {
    pub fn run(self) -> Result<(), Box<dyn std::error::Error>> {
        filter_bedmethyl_by_state(
            &self.input_bed,
            &self.output_bed,
            self.state,
            self.min_gq,
            self.keep_pl,
        )
    }
}

#[derive(Debug, Args)]
pub struct EstimationArgs {
    /// Minimum coverage for error-rate estimation
    #[arg(long, default_value_t = 10)]
    pub min_cov: u64,
    /// Maximum coverage for error-rate estimation
    #[arg(long, default_value_t = 100)]
    pub max_cov: u64,
    /// Minimum cov/(cov+d) for error-rate estimation
    #[arg(long, default_value_t = 0.8)]
    pub min_ratio: f64,
}
