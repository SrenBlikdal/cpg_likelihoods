mod common;

use common::{Scratch, write_input};
use cpg_likelihoods::{estimate_error_c2n_from_bedmethyl, likelihoods_with_normalization};

#[test]
fn counts_are_checked_and_equal_coverage_bounds_are_allowed() {
    let scratch = Scratch::new();
    let input = scratch.path("input.bed");
    write_input(&input);
    let e = estimate_error_c2n_from_bedmethyl(&input, 30, 30, 0.8).unwrap();
    assert!((e - 2.0 / 32.0).abs() < 1e-12);
}

#[test]
fn all_readme_example_rows_match_the_model() {
    for (cov, d, expected) in [
        (30, 2, [0, 62, 409]),
        (28, 0, [0, 80, 409]),
        (15, 16, [145, 0, 131]),
        (13, 16, [151, 0, 107]),
        (1, 29, [409, 71, 0]),
        (30, 8, [8, 0, 330]),
    ] {
        let result = likelihoods_with_normalization(cov, d, 0.03333, 2.0).unwrap();
        assert_eq!(
            [result.pl_hom_cpg, result.pl_het_cpg, result.pl_non_cpg],
            expected
        );
    }
}
