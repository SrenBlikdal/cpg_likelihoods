mod cli;
mod summarize_likelihoods;

use clap::Parser;
use cli::Cli;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    Cli::parse().command.run()
}
