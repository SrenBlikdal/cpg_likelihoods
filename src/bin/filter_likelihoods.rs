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
            .replace('.', "")
            .to_ascii_lowercase();

        match normalized.as_str() {
            "cpg_hom" | "hom_cpg" | "hom" => Ok(Self::HomCpg),
            "cpg_het" | "het_cpg" | "het" => Ok(Self::HetCpg),
            "non_cpg" | "noncpg" | "non" => Ok(Self::NonCpg),
            "cpg_hom_het" | "cpg_homhet" | "hom_het_cpg" | "homhet" | "hom_het"
            => Ok(Self::HomHetCpg),
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

pub fn filter_bedmethyl_by_state<P: AsRef<Path>>(
    input_bed: P,
    output_bed: P,
    state: FilterState,
    min_gq: u32,
    keep_pl: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let reader = open_bedmethyl_reader(input_bed)?;
    let mut writer = BedWriter::create(output_bed)?;

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            writeln!(writer)?;
            continue;
        }

        if line.starts_with('#') {
            writeln!(writer, "{line}")?;
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 21 {
            continue;
        }

        let pl_hom: u32 = match fields[18].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let pl_het: u32 = match fields[19].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let pl_non: u32 = match fields[20].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let locus = classify_from_pl(pl_hom, pl_het, pl_non, min_gq);
        if !state.matches_class(locus) {
            continue;
        }

        if keep_pl {
            writeln!(writer, "{line}")?;
        } else {
            let filtered_fields = &fields[..fields.len().saturating_sub(3)];
            writeln!(writer, "{}", filtered_fields.join("\t"))?;
        }
    }

    writer.finish()
}
