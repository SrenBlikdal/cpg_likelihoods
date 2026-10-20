mod common;

use common::{Scratch, write_input};
use cpg_likelihoods::{
    estimate_error_c2n_from_bedmethyl, likelihoods_with_normalization, parse_site_line,
};

#[test]
fn likelihood_parameters_must_be_finite_and_in_range() {
    for n in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0, -1.0] {
        assert!(
            likelihoods_with_normalization(30, 2, 0.03, n).is_err(),
            "n={n}"
        );
    }
    for e in [f64::NAN, f64::INFINITY, -0.01, 1.01] {
        assert!(
            likelihoods_with_normalization(30, 2, e, 2.0).is_err(),
            "e={e}"
        );
    }
}

#[test]
fn common_scaling_does_not_change_normalized_likelihoods() {
    let baseline = likelihoods_with_normalization(30, 2, 0.03333, 1.0).unwrap();
    for n in [f64::from_bits(1), 0.5, 2.0, f64::MAX] {
        let result = likelihoods_with_normalization(30, 2, 0.03333, n).unwrap();
        assert_eq!(
            [result.pl_hom_cpg, result.pl_het_cpg, result.pl_non_cpg],
            [
                baseline.pl_hom_cpg,
                baseline.pl_het_cpg,
                baseline.pl_non_cpg
            ]
        );
        assert!((result.post_hom_cpg - baseline.post_hom_cpg).abs() < 1e-12);
    }
    // Rounded values from the README's first worked example.
    assert_eq!(
        [
            baseline.pl_hom_cpg,
            baseline.pl_het_cpg,
            baseline.pl_non_cpg
        ],
        [0, 62, 409]
    );
}

#[test]
fn zero_coverage_and_boundary_error_rates_have_defined_results() {
    for e in [0.0, 0.5, 1.0] {
        let zero = likelihoods_with_normalization(0, 0, e, 2.0).unwrap();
        assert_eq!(
            [zero.pl_hom_cpg, zero.pl_het_cpg, zero.pl_non_cpg],
            [0, 0, 0]
        );
        assert!((zero.post_hom_cpg - 1.0 / 3.0).abs() < 1e-12);
        for (cov, d) in [(30, 0), (0, 30), (30, 2)] {
            let result = likelihoods_with_normalization(cov, d, e, 2.0).unwrap();
            let probabilities = [
                result.post_hom_cpg,
                result.post_het_cpg,
                result.post_non_cpg,
            ];
            assert!(probabilities.iter().all(|p| p.is_finite()));
            assert!((probabilities.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
    }
}

#[test]
fn counts_are_checked_and_equal_coverage_bounds_are_allowed() {
    let scratch = Scratch::new();
    let input = scratch.path("input.bed");
    write_input(&input);
    let e = estimate_error_c2n_from_bedmethyl(&input, 30, 30, 0.8).unwrap();
    assert!((e - 2.0 / 32.0).abs() < 1e-12);
    let data = std::fs::read_to_string(&input).unwrap();
    let mut fields: Vec<_> = data.lines().next().unwrap().split('\t').collect();
    fields[4] = "18446744073709551615";
    fields[17] = "1";
    std::fs::write(&input, fields.join("\t")).unwrap();
    assert!(
        estimate_error_c2n_from_bedmethyl(&input, 1, u64::MAX, 0.0)
            .unwrap()
            .is_finite()
    );
    fields[14] = "18446744073709551615";
    assert!(parse_site_line(&fields.join("\t")).is_err());
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
