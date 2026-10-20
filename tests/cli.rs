mod common;

use common::{Scratch, write_input};
use std::process::{Command, Output};

fn cli(args: &[&std::ffi::OsStr]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_cpg-likelihoods"))
        .args(args)
        .output()
        .unwrap()
}

fn success(output: Output) -> String {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn four_subcommands_complete_a_workflow() {
    let scratch = Scratch::new();
    let input = scratch.path("input with spaces.bed");
    let annotated = scratch.path("annotated.bed");
    let filtered = scratch.path("filtered.bed");
    write_input(&input);

    let estimate = success(cli(&["error_c2n".as_ref(), input.as_os_str()]));
    let estimate: f64 = estimate.trim().parse().unwrap();
    assert!((estimate - 2.0 / 32.0).abs() < 1e-9);

    success(cli(&[
        "likelihoods".as_ref(),
        input.as_os_str(),
        "--out".as_ref(),
        annotated.as_os_str(),
    ]));
    let report = success(cli(&[
        "summarize_likelihoods".as_ref(),
        annotated.as_os_str(),
    ]));
    assert!(report.contains("Total loci: 3"));

    success(cli(&[
        "filter_likelihoods".as_ref(),
        "--input-bed".as_ref(),
        annotated.as_os_str(),
        "--output-bed".as_ref(),
        filtered.as_os_str(),
        "--state".as_ref(),
        "CpG_hom/het".as_ref(),
    ]));
    let rows = std::fs::read_to_string(filtered).unwrap();
    assert_eq!(rows.lines().count(), 2);
    assert!(rows.lines().all(|row| row.split('\t').count() == 18));
}

#[test]
fn help_lists_subcommands_and_missing_arguments_fail() {
    let help = success(cli(&["--help".as_ref()]));
    for name in [
        "error_c2n",
        "likelihoods",
        "summarize_likelihoods",
        "filter_likelihoods",
    ] {
        assert!(help.contains(name));
        success(cli(&[name.as_ref(), "--help".as_ref()]));
    }
    assert!(!cli(&["filter_likelihoods".as_ref()]).status.success());
}
