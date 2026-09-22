/// Parse and normalize explicit filter states for likelihood filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterState {
    HomCpg,
    HetCpg,
    NonCpg,
    HomHetCpg,
}

impl FilterState {
    pub fn parse(raw: &str) -> Result<Self, String> {
        let normalized = raw
            .trim()
            .replace(' ', "_")
            .replace('-', "_")
            .replace('/', "_")
            .to_ascii_lowercase();

        match normalized.as_str() {
            "cpg_hom" | "hom_cpg" | "hom" => Ok(Self::HomCpg),
            "cpg_het" | "het_cpg" | "het" => Ok(Self::HetCpg),
            "non_cpg" | "noncpg" | "non_cpg_cpg" => Ok(Self::NonCpg),
            "cpg_hom_het" | "cpg_homhet" | "hom_het_cpg" | "homhet" => Ok(Self::HomHetCpg),
            _ => Err(format!(
                "unsupported state '{raw}'. Use one of: CpG_hom, CpG_het, non-CpG, CpG_hom/het"
            )),
        }
    }

    pub fn matches_class(self, locus: LocusClass) -> bool {
        match self {
            Self::HomCpg => matches!(locus, LocusClass::HomCpg),
            Self::HetCpg => matches!(locus, LocusClass::HetCpg),
            Self::NonCpg => matches!(locus, LocusClass::NonCpg),
            Self::HomHetCpg => matches!(
                locus,
                LocusClass::HomCpg | LocusClass::HetCpg | LocusClass::UnclassifiedCpg
            ),
        }
    }
}

impl std::str::FromStr for FilterState {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}
