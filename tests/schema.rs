mod common;

use common::{Scratch, write_input};
use cpg_likelihoods::{FilterState, append_pl_columns, filter_bedmethyl_by_state};

const HEADER: &str = "#chrom\tstart\tend\tname\tscore\tstrand\tstart2\tend2\tcolor\tN_valid_cov\tpercent\tN_mod\tN_canonical\tN_other_mod\tN_delete\tN_fail\tN_diff\tN_nocall";

#[test]
fn headers_follow_pl_columns_while_comments_and_directives_are_preserved() {
    let scratch = Scratch::new();
    let input = scratch.path("input.bed");
    let annotated = scratch.path("annotated.bed");
    let filtered = scratch.path("filtered.bed");
    write_input(&input);
    let data = std::fs::read_to_string(&input).unwrap();
    let metadata = "# sample comment\ntrack name=sample\nbrowser position chr1:1-10\n";
    std::fs::write(&input, format!("{metadata}{HEADER}\n{data}")).unwrap();
    append_pl_columns(&input, &annotated, 0.03, 2.0).unwrap();
    let result = std::fs::read_to_string(&annotated).unwrap();
    assert!(result.starts_with(metadata));
    let expected_header = format!("{HEADER}\tPL_hom_cpg\tPL_het_cpg\tPL_non_cpg");
    assert_eq!(result.lines().nth(3).unwrap(), expected_header);

    for keep_pl in [false, true] {
        filter_bedmethyl_by_state(&annotated, &filtered, FilterState::HomHetCpg, 20, keep_pl)
            .unwrap();
        let result = std::fs::read_to_string(&filtered).unwrap();
        assert!(result.starts_with(metadata));
        assert_eq!(
            result.lines().nth(3).unwrap(),
            if keep_pl { &expected_header } else { HEADER }
        );
        let expected_columns = if keep_pl { 21 } else { 18 };
        assert!(
            result
                .lines()
                .skip(3)
                .all(|line| line.split('\t').count() == expected_columns)
        );
    }
}

#[test]
fn extra_columns_and_reannotation_fail_instead_of_selecting_wrong_pls() {
    let scratch = Scratch::new();
    let input = scratch.path("input.bed");
    let annotated = scratch.path("annotated.bed");
    let output = scratch.path("output.bed");
    write_input(&input);
    append_pl_columns(&input, &annotated, 0.03, 2.0).unwrap();
    assert!(append_pl_columns(&annotated, &output, 0.03, 2.0).is_err());
    let data = std::fs::read_to_string(&annotated).unwrap();
    let extra = data
        .lines()
        .map(|row| format!("{row}\tunexpected\n"))
        .collect::<String>();
    std::fs::write(&annotated, extra).unwrap();
    assert!(
        filter_bedmethyl_by_state(&annotated, &output, FilterState::HomHetCpg, 20, false).is_err()
    );
}
