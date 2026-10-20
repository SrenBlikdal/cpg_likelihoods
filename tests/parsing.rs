mod common;

use common::{Scratch, write_input};
use std::process::Command;

#[test]
fn mixed_delimiters_work_through_every_command() {
    let scratch = Scratch::new();
    let input = scratch.path("mixed.bed");
    let annotated = scratch.path("annotated.bed");
    let filtered = scratch.path("filtered.bed");
    write_input(&input);
    let rows = std::fs::read_to_string(&input).unwrap();
    let mixed = rows
        .lines()
        .map(|row| {
            let fields: Vec<_> = row.split('\t').collect();
            format!("{}\t{}\n", fields[..9].join("\t"), fields[9..].join(" "))
        })
        .collect::<String>();
    std::fs::write(&input, mixed).unwrap();
    for args in [
        vec!["error_c2n".as_ref(), input.as_os_str()],
        vec![
            "likelihoods".as_ref(),
            input.as_os_str(),
            "--out".as_ref(),
            annotated.as_os_str(),
        ],
        vec!["summarize_likelihoods".as_ref(), annotated.as_os_str()],
        vec![
            "filter_likelihoods".as_ref(),
            "--input-bed".as_ref(),
            annotated.as_os_str(),
            "--output-bed".as_ref(),
            filtered.as_os_str(),
            "--state".as_ref(),
            "homhet".as_ref(),
        ],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_cpg-likelihoods"))
            .args(args)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let output = std::fs::read_to_string(annotated).unwrap();
    assert!(output.lines().all(|line| line.split('\t').count() == 21));
    assert_eq!(
        std::fs::read_to_string(filtered).unwrap().lines().count(),
        2
    );
}
