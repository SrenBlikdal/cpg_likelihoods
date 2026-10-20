mod common;

use common::{Scratch, write_input};
use cpg_likelihoods::append_pl_columns;
use flate2::read::MultiGzDecoder;
use std::io::Read;

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
