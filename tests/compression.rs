mod common;

use common::{Scratch, write_input};
use cpg_likelihoods::{append_pl_columns, open_bedmethyl_reader};
use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use std::io::{Read, Write};

#[test]
fn output_is_really_compressed_and_readable_by_a_gzip_decoder() {
    let scratch = Scratch::new();
    let input = scratch.path("input.bed");
    write_input(&input);
    for extension in ["gz", "bgz", "bgzf"] {
        let output = scratch.path(&format!("output.bed.{extension}"));
        append_pl_columns(&input, &output, 0.03, 2.0).unwrap();
        let bytes = std::fs::read(output).unwrap();
        assert_eq!(&bytes[..2], &[0x1f, 0x8b], "{extension}");
        if extension != "gz" {
            assert_eq!(&bytes[12..14], b"BC", "BGZF extra header");
        }
        let mut text = String::new();
        MultiGzDecoder::new(bytes.as_slice())
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text.lines().count(), 3);
        assert!(text.lines().all(|line| line.split('\t').count() == 21));
    }
}

#[test]
fn reader_accepts_independent_gzip_bgzf_and_concatenated_members() {
    let scratch = Scratch::new();
    let input = scratch.path("plain.bed");
    write_input(&input);
    // More than one BGZF block, and no compression hint in the filename.
    let content = std::fs::read_to_string(&input).unwrap().repeat(2000);
    let gzip = scratch.path("gzip.data");
    let bgzf = scratch.path("bgzf.data");
    let multi = scratch.path("multi.data");

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(content.as_bytes()).unwrap();
    let compressed = encoder.finish().unwrap();
    std::fs::write(&gzip, &compressed).unwrap();
    std::fs::write(
        &multi,
        [compressed.as_slice(), compressed.as_slice()].concat(),
    )
    .unwrap();

    let mut encoder = bgzip::BGZFWriter::new(
        std::fs::File::create(&bgzf).unwrap(),
        bgzip::Compression::default(),
    );
    encoder.write_all(content.as_bytes()).unwrap();
    encoder.close().unwrap();

    for (path, expected) in [
        (&input, std::fs::read_to_string(&input).unwrap()),
        (&gzip, content.clone()),
        (&bgzf, content.clone()),
        (&multi, content.repeat(2)),
    ] {
        let mut actual = String::new();
        open_bedmethyl_reader(path)
            .unwrap()
            .read_to_string(&mut actual)
            .unwrap();
        assert_eq!(actual, expected);
    }
}
