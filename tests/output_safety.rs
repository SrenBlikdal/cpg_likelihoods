mod common;

use common::{Scratch, write_input};
use cpg_likelihoods::{FilterState, append_pl_columns, filter_bedmethyl_by_state};

#[test]
fn annotation_rejects_same_file_and_hard_link_without_changing_input() {
    let scratch = Scratch::new();
    let input = scratch.path("input.bed");
    let alias = scratch.path("hard-link.bed");
    write_input(&input);
    std::fs::hard_link(&input, &alias).unwrap();
    let original = std::fs::read(&input).unwrap();
    for output in [&input, &alias] {
        assert!(append_pl_columns(&input, output, 0.03, 2.0).is_err());
        assert_eq!(std::fs::read(&input).unwrap(), original);
    }
}

#[test]
fn filtering_rejects_same_file_without_changing_input() {
    let scratch = Scratch::new();
    let input = scratch.path("input.bed");
    write_input(&input);
    let original = std::fs::read(&input).unwrap();
    assert!(filter_bedmethyl_by_state(&input, &input, FilterState::HomCpg, 20, false).is_err());
    assert_eq!(std::fs::read(&input).unwrap(), original);
}
