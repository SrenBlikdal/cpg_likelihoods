use crate::Site;

pub const BED_COLUMNS: usize = 18;
pub const ANNOTATED_COLUMNS: usize = 21;

fn header_fields(line: &str) -> Option<(bool, Vec<&str>)> {
    let line = line.trim();
    let commented = line.starts_with('#');
    let fields: Vec<_> = line
        .strip_prefix('#')
        .unwrap_or(line)
        .split_ascii_whitespace()
        .collect();
    if fields.first() == Some(&"chrom") && matches!(fields.get(1), Some(&"start" | &"chromStart")) {
        Some((commented, fields))
    } else {
        None
    }
}

/// Borrow fields and distinguish metadata from malformed records.
pub struct BedMethylRecord<'a> {
    pub fields: Vec<&'a str>,
}

impl<'a> BedMethylRecord<'a> {
    pub fn parse(line: &'a str) -> Result<Option<Self>, String> {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track ")
            || line.starts_with("browser ")
        {
            return Ok(None);
        }
        let fields: Vec<_> = line.split_ascii_whitespace().collect();
        if header_fields(line).is_some() {
            return Ok(None);
        }
        if fields.len() < BED_COLUMNS {
            return Err(format!(
                "expected at least {BED_COLUMNS} bedMethyl columns, found {}",
                fields.len()
            ));
        }
        Ok(Some(Self { fields }))
    }

    pub fn site(&self) -> Result<Site, String> {
        let fields = &self.fields;
        let parse_count = |index: usize, name: &str| -> Result<u64, String> {
            fields
                .get(index)
                .ok_or_else(|| format!("missing {name}"))?
                .parse()
                .map_err(|_| {
                    format!(
                        "column {} ({name}): expected an unsigned integer, found '{}'",
                        index + 1,
                        fields[index]
                    )
                })
        };
        let cov = parse_count(4, "coverage")?;
        let deleted = parse_count(14, "N_delete")?;
        let different = parse_count(16, "N_diff")?;
        let no_call = parse_count(17, "N_nocall")?;
        let d = deleted + different + no_call;
        Ok(Site { cov, d })
    }

    pub fn pl(&self) -> Result<[u32; 3], String> {
        if self.fields.len() < ANNOTATED_COLUMNS {
            return Err(format!(
                "expected at least {ANNOTATED_COLUMNS} PL-annotated columns, found {}",
                self.fields.len()
            ));
        }
        let mut pl = [0; 3];
        for (i, value) in pl.iter_mut().enumerate() {
            *value = self.fields[18 + i].parse().map_err(|_| {
                format!(
                    "column {} (PL): expected an unsigned 32-bit integer, found '{}'",
                    19 + i,
                    self.fields[18 + i]
                )
            })?;
        }
        Ok(pl)
    }
}
