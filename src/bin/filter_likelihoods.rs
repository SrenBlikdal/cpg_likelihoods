use clap::Parser;
use cpg_likelihoods::{filter_bedmethyl_by_state, FilterState};

#[derive(Debug, Parser)]
#[command(name = "filter_likelihoods")]
#[command(about = "Filter PL-annotated bedMethyl rows by CpG state")]
struct Cli {
    #[arg(long)]
    input_bed: String,

    #[arg(long)]
    output_bed: String,

    #[arg(long)]
    state: FilterState,

    #[arg(long, default_value_t = 20)]
    min_gq: u32,

    #[arg(long, default_value_t = false)]
    keep_pl: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    filter_bedmethyl_by_state(
        &cli.input_bed,
        &cli.output_bed,
        cli.state,
        cli.min_gq,
        cli.keep_pl,
    )
}
