mod common;

use common::{Scratch, write_input};
use std::process::Command;

#[test]
fn count_errors_name_the_field_and_metadata_is_distinct_from_invalid_data() {
    for metadata in [
        "",
        "  ",
        "# comment",
        "  # comment",
        "track name=sample",
        "browser position chr1:1-10",
    ] {
        assert!(
            cpg_likelihoods::parse_site_line(metadata)
                .unwrap()
                .is_none()
        );
    }
    let scratch = Scratch::new();
    let input = scratch.path("input.bed");
    write_input(&input);
    let data = std::fs::read_to_string(input).unwrap();
    let mut fields: Vec<_> = data.lines().next().unwrap().split('\t').collect();
    fields[17] = "not-a-count";
    let error = cpg_likelihoods::parse_site_line(&fields.join(" ")).unwrap_err();
    assert!(error.contains("column 18 (N_nocall)"), "{error}");
}

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

#[test]
fn malformed_data_fails_with_path_and_line_number_in_all_commands() {
    let scratch = Scratch::new();
    let input = scratch.path("malformed.bed");
    let out = scratch.path("output.bed");
    std::fs::write(&input, "# a comment\nchr1\tbroken\n").unwrap();
    for args in [
        vec!["error_c2n".as_ref(), input.as_os_str()],
        vec![
            "likelihoods".as_ref(),
            input.as_os_str(),
            "--out".as_ref(),
            out.as_os_str(),
            "--e-c2n".as_ref(),
            "0.03".as_ref(),
        ],
        vec!["summarize_likelihoods".as_ref(), input.as_os_str()],
        vec![
            "filter_likelihoods".as_ref(),
            "--input-bed".as_ref(),
            input.as_os_str(),
            "--output-bed".as_ref(),
            out.as_os_str(),
            "--state".as_ref(),
            "homhet".as_ref(),
        ],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_cpg-likelihoods"))
            .args(args)
            .output()
            .unwrap();
        assert!(!output.status.success());
        let error = String::from_utf8_lossy(&output.stderr);
        assert!(error.contains("malformed.bed:2:"), "{error}");
        assert!(error.contains("18"), "{error}");
    }
    assert!(
        !out.exists(),
        "failed processing must not publish a partial output"
    );
}
