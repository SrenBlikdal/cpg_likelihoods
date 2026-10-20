use crate::Site;
use std::io::{self, Write};

pub const BED_COLUMNS: usize = 18;
pub const ANNOTATED_COLUMNS: usize = 21;
const PL_NAMES: [&str; 3] = ["PL_hom_cpg", "PL_het_cpg", "PL_non_cpg"];

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

fn require_columns(actual: usize, expected: usize) -> Result<(), String> {
    if actual != expected {
        Err(format!(
            "expected exactly {expected} columns, found {actual}"
        ))
    } else {
        Ok(())
    }
}

fn write_header<'a>(
    writer: &mut impl Write,
    commented: bool,
    fields: impl IntoIterator<Item = &'a str>,
) -> io::Result<()> {
    if commented {
        writer.write_all(b"#")?;
    }
    for (index, field) in fields.into_iter().enumerate() {
        if index > 0 {
            writer.write_all(b"\t")?;
        }
        writer.write_all(field.as_bytes())?;
    }
    writer.write_all(b"\n")
}

pub(crate) fn append_pl_header(writer: &mut impl Write, line: &str) -> io::Result<()> {
    let Some((commented, fields)) = header_fields(line) else {
        return writeln!(writer, "{line}");
    };
    require_columns(fields.len(), BED_COLUMNS)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    write_header(writer, commented, fields.into_iter().chain(PL_NAMES))
}

pub(crate) fn filter_pl_header(
    writer: &mut impl Write,
    line: &str,
    keep_pl: bool,
) -> io::Result<()> {
    let Some((commented, fields)) = header_fields(line) else {
        return writeln!(writer, "{line}");
    };
    require_columns(fields.len(), ANNOTATED_COLUMNS)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    if fields[BED_COLUMNS..] != PL_NAMES {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "expected PL_hom_cpg, PL_het_cpg, PL_non_cpg as header columns 19-21",
        ));
    }
    let columns = if keep_pl {
        ANNOTATED_COLUMNS
    } else {
        BED_COLUMNS
    };
    write_header(writer, commented, fields[..columns].iter().copied())
}

/// A modkit bedMethyl row. Writers normalize the parsed fields to tabs.
pub struct BedMethylRecord<'a> {
    pub fields: Vec<&'a str>,
    pub site: Site,
}

impl<'a> BedMethylRecord<'a> {
    /// Blank lines, comments, browser directives, and column headers return None.
    /// Invalid data returns an error rather than looking like a skipped comment.
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
        let d = deleted
            .checked_add(different)
            .and_then(|d| d.checked_add(no_call))
            .ok_or_else(|| "non-CpG counts exceed u64 range".to_owned())?;
        Ok(Some(Self {
            fields,
            site: Site { cov, d },
        }))
    }

    pub fn pl(&self) -> Result<[u32; 3], String> {
        self.require_columns(ANNOTATED_COLUMNS)?;
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

    pub fn require_columns(&self, expected: usize) -> Result<(), String> {
        require_columns(self.fields.len(), expected)
    }
}
