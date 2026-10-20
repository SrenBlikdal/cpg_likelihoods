use crate::Site;

pub const BED_COLUMNS: usize = 18;
pub const ANNOTATED_COLUMNS: usize = 21;

/// Borrow fields once and share the existing count and PL parsers.
pub struct BedMethylRecord<'a> {
    pub fields: Vec<&'a str>,
}

impl<'a> BedMethylRecord<'a> {
    pub fn parse(line: &'a str) -> Option<Self> {
        if line.is_empty() || line.starts_with('#') {
            return None;
        }
        Some(Self {
            fields: line.split_ascii_whitespace().collect(),
        })
    }

    pub fn site(&self) -> Option<Site> {
        if self.fields.len() < BED_COLUMNS {
            return None;
        }
        let cov = self.fields.get(4)?.parse().ok()?;
        let deleted: u64 = self.fields.get(14)?.parse().ok()?;
        let different: u64 = self.fields.get(16)?.parse().ok()?;
        let no_call: u64 = self.fields.get(17)?.parse().ok()?;
        Some(Site {
            cov,
            d: deleted + different + no_call,
        })
    }

    pub fn pl(&self) -> Option<[u32; 3]> {
        if self.fields.len() < ANNOTATED_COLUMNS {
            return None;
        }
        Some([
            self.fields[18].parse().ok()?,
            self.fields[19].parse().ok()?,
            self.fields[20].parse().ok()?,
        ])
    }
}
