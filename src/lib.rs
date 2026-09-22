use bgzip::{BGZFReader, BGZFWriter, Compression as BgzipCompression};
use flate2::Compression as GzipCompression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BedCompression {
    Plain,
    Gzip,
    Bgzip,
}

fn detect_bed_compression(path: &Path) -> BedCompression {
    let path_string = path.to_string_lossy();
    let path_lower = path_string.to_ascii_lowercase();
    if path_lower.ends_with(".bgz") {
        BedCompression::Bgzip
    } else if path_lower.ends_with(".gz") {
        BedCompression::Gzip
    } else {
        BedCompression::Plain
    }
}

pub fn open_bedmethyl_reader<P: AsRef<Path>>(
    path: P,
) -> Result<Box<dyn BufRead>, Box<dyn std::error::Error>> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref)?;

    match detect_bed_compression(path_ref) {
        BedCompression::Plain => Ok(Box::new(BufReader::new(file))),
        BedCompression::Gzip => Ok(Box::new(BufReader::new(MultiGzDecoder::new(file)))),
        BedCompression::Bgzip => Ok(Box::new(BGZFReader::new(file)?)),
    }
}

enum BedWriter {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<File>),
    Bgzip(BGZFWriter<File>),
}

impl BedWriter {
    fn create<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let path_ref = path.as_ref();
        let file = File::create(path_ref)?;

        let writer = match detect_bed_compression(path_ref) {
            BedCompression::Plain => Self::Plain(BufWriter::new(file)),
            BedCompression::Gzip => Self::Gzip(GzEncoder::new(file, GzipCompression::default())),
            BedCompression::Bgzip => Self::Bgzip(BGZFWriter::new(file, BgzipCompression::default())),
        };
        Ok(writer)
    }

    fn finish(self) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            Self::Plain(mut writer) => {
                writer.flush()?;
                Ok(())
            }
            Self::Gzip(writer) => {
                let mut file = writer.finish()?;
                file.flush()?;
                Ok(())
            }
            Self::Bgzip(writer) => {
                writer.close()?;
                Ok(())
            }
        }
    }
}

impl Write for BedWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buf),
            Self::Gzip(writer) => writer.write(buf),
            Self::Bgzip(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
            Self::Bgzip(writer) => writer.flush(),
        }
    }
}

/// One parsed bedMethyl locus using your fixed schema:
/// - cov = col5 (1-based)
/// - d   = col15 + col17 + col18 (1-based)
#[derive(Debug, Clone, Copy)]
pub struct Site {
    pub cov: u64,
    pub d: u64,
}

/// Output bundle from likelihood computation.
#[derive(Debug, Clone, Copy)]
pub struct LikelihoodResult {
    pub log10_hom_cpg: f64, // g=0
    pub log10_het_cpg: f64, // g=1
    pub log10_non_cpg: f64, // g=2
    pub post_hom_cpg: f64,  // normalized probs
    pub post_het_cpg: f64,
    pub post_non_cpg: f64,
    pub pl_hom_cpg: u32,    // GATK-style PL
    pub pl_het_cpg: u32,
    pub pl_non_cpg: u32,
}

/// Parse one data line into Site.
/// Returns None for header/comment/invalid lines.
pub fn parse_site_line(line: &str) -> Option<Site> {
    if line.is_empty() || line.starts_with('#') {
        return None;
    }

    // 0-based indices
    const COV_IDX: usize = 4;   // col5
    const D15_IDX: usize = 14;  // col15
    const D17_IDX: usize = 16;  // col17
    const D18_IDX: usize = 17;  // col18
    const MIN_FIELDS: usize = 18;

    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < MIN_FIELDS {
        return None;
    }

    let cov: u64 = fields[COV_IDX].parse().ok()?;
    let d15: u64 = fields[D15_IDX].parse().ok()?;
    let d17: u64 = fields[D17_IDX].parse().ok()?;
    let d18: u64 = fields[D18_IDX].parse().ok()?;

    Some(Site {
        cov,
        d: d15 + d17 + d18,
    })
}

/// Estimate c2n error rate:
/// ratio = cov/(cov+d)
/// keep: min_cov <= cov <= max_cov and ratio >= min_ratio
/// e_c2n = sum(d_keep) / (sum(d_keep)+sum(cov_keep))
pub fn estimate_error_c2n_from_bedmethyl<P: AsRef<Path>>(
    path: P,
    min_cov: u64,
    max_cov: u64,
    min_ratio: f64,
) -> Result<f64, Box<dyn std::error::Error>> {
    if max_cov <= min_cov {
        return Err("max_cov must be greater than min_cov".into());
    }
    if !(0.0..=1.0).contains(&min_ratio) {
        return Err("min_ratio must be between 0 and 1".into());
    }

    let reader = open_bedmethyl_reader(path)?;

    let mut sum_cov: u128 = 0;
    let mut sum_d: u128 = 0;
    let mut n_sites: u64 = 0;

    for line in reader.lines() {
        let line = line?;
        let Some(site) = parse_site_line(&line) else {
            continue;
        };

        let denom = site.cov + site.d;
        if denom == 0 {
            continue;
        }

        let ratio = site.cov as f64 / denom as f64;
        let keep = site.cov >= min_cov && site.cov <= max_cov && ratio >= min_ratio;
        if keep {
            sum_cov += site.cov as u128;
            sum_d += site.d as u128;
            n_sites += 1;
        }
    }

    if n_sites == 0 {
        return Err("No sites passed filtering criteria".into());
    }

    let total = sum_cov + sum_d;
    if total == 0 {
        return Err("Filtered totals are zero; cannot compute error rate".into());
    }

    Ok(sum_d as f64 / total as f64)
}

fn genotype_log10_likelihood(g: u8, e: f64, cov: u64, d: u64, n: f64) -> Result<f64, &'static str> {
    let (p_mismatch, p_match) = match g {
        0 => (n * e, n * (1.0 - e)), // hom_cpg
        1 => (n / 2.0, n / 2.0),     // het_cpg
        2 => (n * (1.0 - e), n * e), // non_cpg
        _ => return Err("Invalid genotype value. Must be 0,1,2"),
    };

    let mismatch_term = if p_mismatch == 0.0 {
        if d == 0 { 0.0 } else { f64::NEG_INFINITY }
    } else {
        (d as f64) * p_mismatch.log10()
    };

    let match_term = if p_match == 0.0 {
        if cov == 0 { 0.0 } else { f64::NEG_INFINITY }
    } else {
        (cov as f64) * p_match.log10()
    };

    Ok(mismatch_term + match_term)
}

fn log10_to_pl(log10s: [f64; 3]) -> [u32; 3] {
    let max_l = log10s.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let mut out = [0_u32; 3];
    for (i, &li) in log10s.iter().enumerate() {
        let delta = max_l - li; // >=0
        if !delta.is_finite() {
            out[i] = u32::MAX;
            continue;
        }
        let v = (10.0 * delta).round();
        out[i] = if v < 0.0 {
            0
        } else if v > u32::MAX as f64 {
            u32::MAX
        } else {
            v as u32
        };
    }
    out
}

fn normalize_log10(log10s: [f64; 3]) -> [f64; 3] {
    let max_l = log10s.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max_l.is_finite() {
        return [f64::NAN, f64::NAN, f64::NAN];
    }

    let x0 = 10f64.powf(log10s[0] - max_l);
    let x1 = 10f64.powf(log10s[1] - max_l);
    let x2 = 10f64.powf(log10s[2] - max_l);
    let s = x0 + x1 + x2;
    if s == 0.0 || !s.is_finite() {
        return [f64::NAN, f64::NAN, f64::NAN];
    }
    [x0 / s, x1 / s, x2 / s]
}

/// Computes log10 likelihoods + normalized probs + PL values.
pub fn likelihoods_with_normalization(
    cov: u64,
    d: u64,
    e: f64,
    n: f64,
) -> Result<LikelihoodResult, &'static str> {
    if !(0.0..=1.0).contains(&e) {
        return Err("e must be between 0 and 1");
    }
    if n <= 0.0 {
        return Err("n must be > 0");
    }

    let l0 = genotype_log10_likelihood(0, e, cov, d, n)?;
    let l1 = genotype_log10_likelihood(1, e, cov, d, n)?;
    let l2 = genotype_log10_likelihood(2, e, cov, d, n)?;

    let log10s = [l0, l1, l2];
    let post = normalize_log10(log10s);
    let pl = log10_to_pl(log10s);

    Ok(LikelihoodResult {
        log10_hom_cpg: l0,
        log10_het_cpg: l1,
        log10_non_cpg: l2,
        post_hom_cpg: post[0],
        post_het_cpg: post[1],
        post_non_cpg: post[2],
        pl_hom_cpg: pl[0],
        pl_het_cpg: pl[1],
        pl_non_cpg: pl[2],
    })
}

/// Read input bed and append PL columns as new cols 19-21 (1-based).
pub fn append_pl_columns<P: AsRef<Path>>(
    input_bed: P,
    output_bed: P,
    e_c2n: f64,
    n: f64,
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
            writeln!(writer, "{}\tPL_hom_cpg\tPL_het_cpg\tPL_non_cpg", line)?;
            continue;
        }

        if let Some(site) = parse_site_line(&line) {
            let lk = likelihoods_with_normalization(site.cov, site.d, e_c2n, n)
                .map_err(|e| format!("likelihood error: {e}"))?;
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                line, lk.pl_hom_cpg, lk.pl_het_cpg, lk.pl_non_cpg
            )?;
        } else {
            // header-like/non-data row without '#'
            writeln!(writer, "{}\tPL_hom_cpg\tPL_het_cpg\tPL_non_cpg", line)?;
        }
    }

    writer.finish()
}

/// Classification result for summary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocusClass {
    HomCpg,
    HetCpg,
    NonCpg,
    UnclassifiedCpg,
    Unclassified,
}

/// Choose best class from PLs with threshold rules:
/// - Unique minimum PL required.
/// - GQ = second_smallest - smallest must be >= min_gq.
pub fn classify_from_pl(
    pl_hom: u32,
    pl_het: u32,
    pl_non: u32,
    min_gq: u32,
) -> LocusClass {
    let mut vals = [(0_u8, pl_hom), (1_u8, pl_het), (2_u8, pl_non)];
    vals.sort_by_key(|x| x.1);

    let best = vals[0];
    let second = vals[1];

    // Tie at best PL -> unclassified
    if best.1 == second.1 {
        return LocusClass::Unclassified;
    }

    let gq = second.1 - best.1;

    // CpG-rescue bucket: best is hom/het, low GQ overall,
    // but non_cpg is still clearly worse than best CpG.
    if (best.0 == 0 || best.0 == 1) && gq < min_gq {
        let best_cpg_pl = pl_hom.min(pl_het);
        let cpg_vs_non_delta = pl_non.saturating_sub(best_cpg_pl);
        if cpg_vs_non_delta >= min_gq {
            return LocusClass::UnclassifiedCpg;
        }
        return LocusClass::Unclassified;
    }

    if gq < min_gq {
        return LocusClass::Unclassified;
    }

    match best.0 {
        0 => LocusClass::HomCpg,
        1 => LocusClass::HetCpg,
        2 => LocusClass::NonCpg,
        _ => LocusClass::Unclassified,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::{SystemTime, UNIX_EPOCH};

    static TEST_COUNTER: AtomicU64 = AtomicU64::new(0);

    struct Cleanup {
        paths: Vec<std::path::PathBuf>,
    }

    impl Cleanup {
        fn new(paths: Vec<std::path::PathBuf>) -> Self {
            Self { paths }
        }
    }

    impl Drop for Cleanup {
        fn drop(&mut self) {
            for path in &self.paths {
                let _ = fs::remove_file(path);
            }
        }
    }

    fn temp_path(filename: &str) -> std::path::PathBuf {
        let suffix = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock should be after UNIX_EPOCH")
            .as_nanos();
        let counter = TEST_COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!("cpg_likelihoods_{suffix}_{counter}_{filename}"))
    }

    fn sample_bed_lines() -> Vec<&'static str> {
        vec![
            "#chrom\tstart\tend\tname\tcov\tscore\tstrand\ta\tb\tc\td\te\tf\tg\td15\ti\td17\td18",
            "chr1\t0\t1\tlocus1\t10\t0\t+\t0\t0\t0\t0\t0\t0\t0\t2\t0\t1\t1",
        ]
    }

    fn write_lines(path: &Path, lines: &[&str]) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = BedWriter::create(path)?;
        for line in lines {
            writeln!(writer, "{line}")?;
        }
        writer.finish()
    }

    #[test]
    fn estimate_error_c2n_reads_plain_gz_and_bgz() -> Result<(), Box<dyn std::error::Error>> {
        let lines = sample_bed_lines();
        let plain_path = temp_path("input.bed");
        let gzip_path = temp_path("input.bed.gz");
        let bgzip_path = temp_path("input.bed.bgz");
        let _cleanup = Cleanup::new(vec![plain_path.clone(), gzip_path.clone(), bgzip_path.clone()]);

        write_lines(&plain_path, &lines)?;
        write_lines(&gzip_path, &lines)?;
        write_lines(&bgzip_path, &lines)?;

        let expected = 4.0 / 14.0;
        let plain = estimate_error_c2n_from_bedmethyl(&plain_path, 1, 100, 0.0)?;
        let gz = estimate_error_c2n_from_bedmethyl(&gzip_path, 1, 100, 0.0)?;
        let bgz = estimate_error_c2n_from_bedmethyl(&bgzip_path, 1, 100, 0.0)?;

        assert!((plain - expected).abs() < 1e-12);
        assert!((gz - expected).abs() < 1e-12);
        assert!((bgz - expected).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn append_pl_columns_writes_gz_and_bgz() -> Result<(), Box<dyn std::error::Error>> {
        let lines = sample_bed_lines();
        let input_path = temp_path("append_input.bed");
        let output_gz_path = temp_path("append_output.bed.gz");
        let output_bgz_path = temp_path("append_output.bed.bgz");
        let _cleanup = Cleanup::new(vec![
            input_path.clone(),
            output_gz_path.clone(),
            output_bgz_path.clone(),
        ]);

        write_lines(&input_path, &lines)?;
        append_pl_columns(&input_path, &output_gz_path, 0.1, 2.0)?;
        append_pl_columns(&input_path, &output_bgz_path, 0.1, 2.0)?;

        let gz_lines: Vec<String> = open_bedmethyl_reader(&output_gz_path)?.lines().collect::<Result<_, _>>()?;
        let bgz_lines: Vec<String> = open_bedmethyl_reader(&output_bgz_path)?.lines().collect::<Result<_, _>>()?;

        for out_lines in [&gz_lines, &bgz_lines] {
            assert_eq!(out_lines.len(), 2);
            assert!(out_lines[0].ends_with("\tPL_hom_cpg\tPL_het_cpg\tPL_non_cpg"));
            assert_eq!(out_lines[1].split('\t').count(), 21);
        }

        Ok(())
    }
}

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