use assert_cmd::Command;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::Path;
use tempfile::TempDir;

fn random_seq(len: usize, rng: &mut StdRng) -> String {
    (0..len)
        .map(|_| b"ACGT"[rng.gen_range(0..4)] as char)
        .collect()
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            _ => 'C',
        })
        .collect()
}

/// (name, sequence, quality character) records as FASTQ text.
fn fastq(records: &[(&str, &str, char)]) -> String {
    records
        .iter()
        .map(|(name, seq, q)| format!("@{name}\n{seq}\n+\n{}\n", q.to_string().repeat(seq.len())))
        .collect()
}

fn read_maybe_gz(path: &Path) -> String {
    let bytes = fs::read(path).unwrap();
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut text = String::new();
        MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .unwrap();
        text
    } else {
        String::from_utf8(bytes).unwrap()
    }
}

#[test]
fn deduplicates_gzipped_fastq_and_writes_assignments() {
    let mut rng = StdRng::seed_from_u64(1);
    let (a, b) = (random_seq(1500, &mut rng), random_seq(1200, &mut rng));
    let a_rc = reverse_complement(&a);
    let input = fastq(&[
        ("r0", &a, '5'),
        ("r1", &b, '5'),
        ("r2", &a_rc, '?'), // best copy of `a`, opposite strand
        ("r3", &a, '+'),
    ]);
    let dir = TempDir::new().unwrap();
    let (in_path, out_path, assign_path) = (
        dir.path().join("in.fastq.gz"),
        dir.path().join("out.fastq.gz"),
        dir.path().join("assign.tsv"),
    );
    let mut gz = GzEncoder::new(File::create(&in_path).unwrap(), Compression::default());
    gz.write_all(input.as_bytes()).unwrap();
    gz.finish().unwrap();

    Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg(&in_path)
        .arg(&out_path)
        .arg("--assignments")
        .arg(&assign_path)
        .assert()
        .success();

    let output = read_maybe_gz(&out_path);
    let names: Vec<_> = output.lines().step_by(4).collect();
    assert_eq!(names, ["@r1", "@r2"]);
    assert_eq!(
        fs::read_to_string(&assign_path).unwrap(),
        "read_idx\texemplar_idx\n0\t2\n1\t1\n2\t2\n3\t2\n"
    );
}

#[test]
fn no_rc_keeps_opposite_strands_apart_and_plain_output() {
    let mut rng = StdRng::seed_from_u64(2);
    let a = random_seq(1000, &mut rng);
    let input = fastq(&[("r0", &a, '5'), ("r1", &reverse_complement(&a), '5')]);
    let dir = TempDir::new().unwrap();
    let (in_path, out_path) = (dir.path().join("in.fastq"), dir.path().join("out.fastq"));
    fs::write(&in_path, input).unwrap();

    Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg(&in_path)
        .arg(&out_path)
        .arg("--no-rc")
        .assert()
        .success();
    assert_eq!(read_maybe_gz(&out_path).lines().count(), 8);
}

fn run_on(input: &str, extra: &[&str]) -> (assert_cmd::assert::Assert, TempDir) {
    let dir = TempDir::new().unwrap();
    let in_path = dir.path().join("in.fastq");
    fs::write(&in_path, input).unwrap();
    let assert = Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg(&in_path)
        .arg(dir.path().join("out.fastq"))
        .args(extra)
        .assert();
    (assert, dir)
}

#[test]
fn rejects_malformed_fastq() {
    for input in [
        "@r0\nACGT\n+\nII\n",   // quality length differs
        "r0\nACGT\n+\nIIII\n",  // header without '@'
        "@r0\nACGT\n-\nIIII\n", // separator without '+'
        "@r0\nACGT\n+\n",       // truncated record
        "@r0\nACGT\n",          // truncated record
    ] {
        run_on(input, &[]).0.failure();
    }
}

#[test]
fn accepts_crlf_missing_final_newline_and_empty_input() {
    let (assert, dir) = run_on("@r0 note\r\nACGT\r\n+r0 note\r\nIIII", &[]);
    assert.success();
    // Records are written unchanged apart from line endings.
    assert_eq!(
        fs::read_to_string(dir.path().join("out.fastq")).unwrap(),
        "@r0 note\nACGT\n+r0 note\nIIII\n"
    );
    let (assert, dir) = run_on("", &[]);
    assert.success();
    assert_eq!(
        fs::read_to_string(dir.path().join("out.fastq")).unwrap(),
        ""
    );
}

#[test]
fn gz_output_is_compressed() {
    let mut rng = StdRng::seed_from_u64(3);
    let dir = TempDir::new().unwrap();
    let (in_path, out_path) = (dir.path().join("in.fastq"), dir.path().join("out.fastq.gz"));
    fs::write(&in_path, fastq(&[("r0", &random_seq(500, &mut rng), '5')])).unwrap();
    Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg(&in_path)
        .arg(&out_path)
        .assert()
        .success();
    assert!(fs::read(&out_path).unwrap().starts_with(&[0x1f, 0x8b]));
    assert!(read_maybe_gz(&out_path).starts_with("@r0\n"));
}

#[test]
fn rejects_invalid_error_fraction() {
    for value in ["inf", "NaN", "-0.1", "1.5"] {
        run_on("", &["--max-error-frac", value]).0.failure();
    }
}

#[test]
fn rejects_outputs_that_overwrite_the_input() {
    let dir = TempDir::new().unwrap();
    let in_path = dir.path().join("in.fastq");
    fs::write(&in_path, "@r0\nACGT\n+\nIIII\n").unwrap();
    let out_path = dir.path().join("out.fastq");
    for args in [
        vec![in_path.clone(), in_path.clone()],
        vec![
            in_path.clone(),
            out_path.clone(),
            "--assignments".into(),
            in_path.clone(),
        ],
        vec![
            in_path.clone(),
            out_path.clone(),
            "--assignments".into(),
            out_path.clone(),
        ],
    ] {
        Command::cargo_bin("dedup_ont_fastq")
            .unwrap()
            .args(&args)
            .assert()
            .failure();
    }
    assert_eq!(
        fs::read_to_string(&in_path).unwrap(),
        "@r0\nACGT\n+\nIIII\n"
    );
}
