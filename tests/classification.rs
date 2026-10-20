use cpg_likelihoods::{FilterState, LocusClass, classify_from_pl};

#[test]
fn a_cpg_tie_can_be_confidently_classified_as_combined_cpg() {
    let class = classify_from_pl(0, 0, 100, 20);
    assert_eq!(class, LocusClass::UnclassifiedCpg);
    assert!(FilterState::HomHetCpg.matches_class(class));
    assert!(!FilterState::HomCpg.matches_class(class));
    assert!(!FilterState::HetCpg.matches_class(class));
    assert_eq!(classify_from_pl(0, 0, 100, 0), class);
}

#[test]
fn ties_involving_non_cpg_remain_unclassified_even_at_zero_threshold() {
    for pl in [[0, 0, 0], [0, 100, 0], [100, 0, 0]] {
        for threshold in [0, 20] {
            assert_eq!(
                classify_from_pl(pl[0], pl[1], pl[2], threshold),
                LocusClass::Unclassified
            );
        }
    }
}

#[test]
fn classification_obeys_quality_boundaries_for_separate_and_combined_states() {
    assert_eq!(classify_from_pl(0, 20, 100, 20), LocusClass::HomCpg);
    assert_eq!(classify_from_pl(20, 0, 100, 20), LocusClass::HetCpg);
    assert_eq!(classify_from_pl(100, 20, 0, 20), LocusClass::NonCpg);
    assert_eq!(classify_from_pl(0, 19, 20, 20), LocusClass::UnclassifiedCpg);
    assert_eq!(classify_from_pl(0, 0, 19, 20), LocusClass::Unclassified);
    assert_eq!(classify_from_pl(0, 0, 20, 20), LocusClass::UnclassifiedCpg);
}
