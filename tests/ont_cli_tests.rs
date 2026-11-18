use assert_cmd::{Command, assert::Assert};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::borrow::Cow;
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
    String::from_utf8(ont_sim::reverse_complement(seq.as_bytes())).unwrap()
}

fn stderr_of(result: &Assert) -> Cow<'_, str> {
    String::from_utf8_lossy(&result.get_output().stderr)
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut gz = GzEncoder::new(Vec::new(), Compression::default());
    gz.write_all(bytes).unwrap();
    gz.finish().unwrap()
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
    fs::write(&in_path, gzip(input.as_bytes())).unwrap();

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
        .current_dir(dir.path())
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
fn cli_exposes_user_decisions() {
    let help = Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg("--help")
        .assert()
        .success();
    let text = String::from_utf8_lossy(&help.get_output().stdout);
    for flag in [
        "--assignments",
        "--protect-divergence",
        "--noise-scale",
        "--no-rc",
        "--max-end-offset",
        "--assumed-error-rate",
    ] {
        assert!(text.contains(flag), "{text}");
    }
    assert!(text.contains("regular file"), "{text}");
    let noise_help = text.split_once("--noise-scale").unwrap().1;
    let noise_help = noise_help.split("\n      --").next().unwrap();
    assert!(noise_help.contains("[default: auto]"), "{text}");
    // Tuning constants stay library-only; the CLI exposes user decisions.
    for removed in [
        "--end-check-len",
        "--max-end-errors",
        "--max-error-frac",
        "--quality-calibration",
        "--quality-error-floor",
        "--pair-distance-cv",
        "--seed-drift-frac",
        "--kmer-len",
        "--window-len",
        "--num-windows",
    ] {
        assert!(!text.contains(removed), "{removed} in help: {text}");
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

#[test]
fn plain_output_can_go_to_stdout() {
    let mut rng = StdRng::seed_from_u64(4);
    let seq = random_seq(600, &mut rng);
    let dir = TempDir::new().unwrap();
    let in_path = dir.path().join("in.fastq");
    fs::write(&in_path, fastq(&[("r0", &seq, '5'), ("r1", &seq, '5')])).unwrap();
    let output = Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg(&in_path)
        .arg("/dev/stdout")
        .output()
        .unwrap();
    assert!(output.status.success());
    assert_eq!(String::from_utf8(output.stdout).unwrap().lines().count(), 4);
}

#[test]
fn rejects_output_hard_linked_to_input() {
    let dir = TempDir::new().unwrap();
    let in_path = dir.path().join("in.fastq");
    fs::write(&in_path, "@r0\nACGT\n+\nIIII\n").unwrap();
    let link = dir.path().join("link.fastq");
    fs::hard_link(&in_path, &link).unwrap();
    Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg(&in_path)
        .arg(&link)
        .assert()
        .failure();
    assert_eq!(
        fs::read_to_string(&in_path).unwrap(),
        "@r0\nACGT\n+\nIIII\n"
    );
}

#[test]
fn detects_gzip_by_content_and_respects_endpoint_tolerance() {
    // Two reads whose starts differ by 20 bases: separate by default,
    // merged once the endpoint tolerance covers the shift. The input is
    // gzip-compressed despite its plain .fastq name.
    let mut rng = StdRng::seed_from_u64(5);
    let genome = random_seq(2200, &mut rng);
    let input = fastq(&[
        ("r0", &genome[100..2100], '5'),
        ("r1", &genome[120..2100], '5'),
    ]);
    let dir = TempDir::new().unwrap();
    let in_path = dir.path().join("in.fastq");
    fs::write(&in_path, gzip(input.as_bytes())).unwrap();
    let count = |extra: &[&str]| {
        let out = dir.path().join("out.fastq");
        Command::cargo_bin("dedup_ont_fastq")
            .unwrap()
            .arg(&in_path)
            .arg(&out)
            .args(extra)
            .assert()
            .success();
        read_maybe_gz(&out).lines().count() / 4
    };
    assert_eq!(count(&[]), 2);
    assert_eq!(count(&["--max-end-offset", "25"]), 1);
}

#[test]
fn names_bad_records() {
    let (assert, _dir) = run_on("@r0\nACGT\n+\nIIII\n@r1\nACGT\n+\nII\n", &[]);
    let stderr = stderr_of(&assert.failure()).into_owned();
    assert!(stderr.contains("FASTQ record 2"), "{stderr}");
    assert!(stderr.contains("in.fastq"), "{stderr}");
}

#[test]
fn rejects_aliased_new_outputs_and_shared_streams() {
    let dir = TempDir::new().unwrap();
    let in_path = dir.path().join("in.fastq");
    fs::write(&in_path, "@r0\nACGT\n+\nIIII\n").unwrap();
    fs::create_dir(dir.path().join("sub")).unwrap();
    let cases = [
        // Same new file reached through `sub/..`.
        (
            dir.path().join("out.fastq"),
            dir.path().join("sub/../out.fastq"),
        ),
        // Both outputs on one stream.
        ("/dev/stdout".into(), "/proc/self/fd/1".into()),
    ];
    for (out, assignments) in cases {
        Command::cargo_bin("dedup_ont_fastq")
            .unwrap()
            .arg(&in_path)
            .arg(&out)
            .arg("--assignments")
            .arg(&assignments)
            .assert()
            .failure();
    }
}

#[test]
fn unreadable_input_error_names_the_path() {
    let dir = TempDir::new().unwrap();
    let assert = Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg(dir.path())
        .arg(dir.path().join("out.fastq"))
        .assert()
        .failure();
    let stderr = stderr_of(&assert);
    assert!(
        stderr.contains(&dir.path().display().to_string()),
        "{stderr}"
    );
}

#[test]
fn validates_statistical_parameters_and_phred_bytes() {
    for (flag, bad, message) in [
        ("--protect-divergence", "NaN", "protect_divergence"),
        ("--protect-divergence", "inf", "protect_divergence"),
        ("--protect-divergence", "0", "protect_divergence"),
        ("--protect-divergence", "1.1", "protect_divergence"),
        ("--assumed-error-rate", "NaN", "assumed_error_rate"),
        ("--assumed-error-rate", "inf", "assumed_error_rate"),
        ("--assumed-error-rate", "1.1", "assumed_error_rate"),
        (
            "--noise-scale",
            "0",
            "expected auto or a finite positive number",
        ),
        (
            "--noise-scale",
            "NaN",
            "expected auto or a finite positive number",
        ),
        (
            "--noise-scale",
            "inf",
            "expected auto or a finite positive number",
        ),
        (
            "--noise-scale",
            "oops",
            "expected auto or a finite positive number",
        ),
    ] {
        run_on("", &[&format!("{flag}={bad}")])
            .0
            .failure()
            .stderr(predicates::str::contains(message));
    }
    for flag in [
        "--protect-divergence",
        "--noise-scale",
        "--assumed-error-rate",
    ] {
        run_on("", &[&format!("{flag}=-0.1")]).0.failure();
    }
    // Absurd multipliers would overflow the noise estimate and silently stop
    // all merging, so they are rejected.
    run_on("", &["--noise-scale", "1000"]).0.failure();
    // Constant invalid qualities are input errors, not unavailable placeholders.
    run_on("@r0\nACGT\n+\n    \n", &[]).0.failure();
    run_on("", &["--assumed-error-rate", "0"]).0.success();
}

#[test]
fn auto_calibration_reports_fallback_and_fixed_scale_bypasses_pilot() {
    let input = "@r0\nACGT\n+\nIJIJ\n";
    let (assert, _) = run_on(input, &[]);
    let stderr = stderr_of(&assert.success()).into_owned();
    for field in [
        "Noise scale: status=fallback(insufficient eligible pairs)",
        "a=1.0000",
        "pairs=0",
        "stop=none",
    ] {
        assert!(stderr.contains(field), "{stderr}");
    }
    for value in ["0.25", "2.5", "6"] {
        let (assert, _) = run_on(input, &["--noise-scale", value]);
        let result = assert.success();
        let stderr = stderr_of(&result);
        assert_eq!(
            stderr.lines().next().unwrap(),
            format!("Noise scale: fixed a={value}")
        );
    }
}

#[test]
fn auto_calibration_changes_duplicate_assignments_on_noisy_reads() {
    use ont_sim::{QualityConfig, add_quality_errors};

    let mut rng = StdRng::seed_from_u64(904);
    // Explicit synthetic miscalibration exercises an upward correction.
    let model = QualityConfig::enabled(0.004, 0.1, 1.0, 2.5, 0.0);
    let mut input = String::new();
    for family in 0..300 {
        let source = random_seq(2400, &mut rng);
        for copy in 0..2 {
            let (sequence, quality, _) =
                add_quality_errors(source.as_bytes(), &model, [0.4, 0.3, 0.3], 1.0, &mut rng);
            input.push_str(&format!(
                "@r{family}_{copy}\n{}\n+\n{}\n",
                String::from_utf8(sequence).unwrap(),
                String::from_utf8(quality).unwrap()
            ));
        }
    }
    let merged_families = |extra: &[&str]| {
        let mut args = vec!["--assignments", "assign.tsv"];
        args.extend_from_slice(extra);
        let (assert, dir) = run_on(&input, &args);
        let result = assert.success();
        let stderr = stderr_of(&result).into_owned();
        let assignments = fs::read_to_string(dir.path().join("assign.tsv")).unwrap();
        let exemplars: Vec<usize> = assignments
            .lines()
            .skip(1)
            .map(|line| line.split_once('\t').unwrap().1.parse().unwrap())
            .collect();
        assert_eq!(exemplars.len(), 600);
        for (idx, &exemplar) in exemplars.iter().enumerate() {
            assert_eq!(idx / 2, exemplar / 2, "unrelated families merged");
        }
        let merged = exemplars.chunks_exact(2).filter(|p| p[0] == p[1]).count();
        (merged, stderr)
    };
    let (automatic, stderr) = merged_families(&[]);
    assert!(stderr.contains("status=estimated"), "{stderr}");
    let (baseline, _) = merged_families(&["--noise-scale", "1"]);
    assert!(
        automatic >= baseline + 30,
        "auto merged {automatic}/300 families, baseline merged {baseline}/300; {stderr}"
    );
    let (fixed, _) = merged_families(&["--noise-scale", "0.5"]);
    assert!(
        automatic >= fixed + 100,
        "auto merged {automatic}/300 families, fixed 0.5 merged {fixed}/300; {stderr}"
    );
}

#[test]
fn pilot_stops_on_target_pairs_replays_and_continues_the_same_reader() {
    let mut rng = StdRng::seed_from_u64(621);
    let mut input = String::new();
    let mut first = String::new();
    let quality = "56".repeat(600);
    for family in 0..500 {
        let sequence = random_seq(1200, &mut rng);
        if family == 0 {
            first = sequence.clone();
        }
        for copy in 0..2 {
            input.push_str(&format!(
                "@r{family}_{copy} note\n{sequence}\n+note\n{quality}\n"
            ));
        }
    }
    let tail = random_seq(1200, &mut rng);
    input.push_str(&format!("@tail_unique\n{tail}\n+\n{quality}\n"));
    let best = format!(
        "@tail_best note\n{first}\n+best note\n{}\n",
        "?@".repeat(600)
    );
    input.push_str(&best);
    let (assert, dir) = run_on(&input, &["--assignments", "assign.tsv"]);
    let stderr = stderr_of(&assert.success()).into_owned();
    for field in [
        "status=estimated",
        "pairs=500; reads=1000",
        "stop=target-pairs",
        "Total reads: 1002",
    ] {
        assert!(stderr.contains(field), "{stderr}");
    }
    let assignments = fs::read_to_string(dir.path().join("assign.tsv")).unwrap();
    let rows: Vec<_> = assignments.lines().skip(1).collect();
    assert_eq!(rows.len(), 1002);
    for (idx, row) in rows.iter().enumerate() {
        let expected = if idx < 2 || idx == 1001 {
            1001
        } else if idx == 1000 {
            1000
        } else {
            idx / 2 * 2
        };
        assert_eq!(*row, format!("{idx}\t{expected}"));
    }
    let output = fs::read_to_string(dir.path().join("out.fastq")).unwrap();
    assert!(output.ends_with(&best));
    assert_eq!(output.lines().count(), 501 * 4);
    let compressed = gzip(input.as_bytes());
    let stdin_assignments = dir.path().join("stdin.tsv");
    let result = Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg("-")
        .arg("--assignments")
        .arg(&stdin_assignments)
        .write_stdin(compressed)
        .assert()
        .success();
    let stdin_stderr = stderr_of(&result);
    assert_eq!(stdin_stderr.lines().next(), stderr.lines().next());
    assert_eq!(fs::read_to_string(stdin_assignments).unwrap(), assignments);
}

#[test]
fn assignments_only_replays_prefix_once_for_files_and_plain_or_gzip_stdin() {
    let mut rng = StdRng::seed_from_u64(811);
    let (a, b) = (random_seq(120, &mut rng), random_seq(120, &mut rng));
    for count in [4, 20_003] {
        let mut input = fastq(&[("first", &a, '5'), ("other", &b, '5')]);
        for idx in 2..count {
            input.push_str(&fastq(&[(
                "copy",
                &a,
                if idx == count - 1 { '?' } else { '5' },
            )]));
        }
        let expected = format!(
            "read_idx\texemplar_idx\n{}",
            (0..count)
                .map(|idx| format!("{idx}\t{}\n", if idx == 1 { 1 } else { count - 1 }))
                .collect::<String>()
        );
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("in.fastq");
        fs::write(&path, &input).unwrap();
        let compressed = gzip(input.as_bytes());
        for fixed in [false, true] {
            for (source, payload) in [
                (path.as_path(), None),
                (Path::new("-"), Some(input.as_bytes())),
                (Path::new("/dev/stdin"), Some(input.as_bytes())),
                (Path::new("-"), Some(compressed.as_slice())),
                (Path::new("/dev/stdin"), Some(compressed.as_slice())),
            ] {
                let assignments = dir.path().join("assign.tsv");
                let mut cmd = Command::cargo_bin("dedup_ont_fastq").unwrap();
                cmd.arg(source).arg("--assignments").arg(&assignments);
                if fixed {
                    cmd.args(["--noise-scale", "1"]);
                }
                if let Some(payload) = payload {
                    cmd.write_stdin(payload);
                }
                let result = cmd.assert().success();
                let stderr = stderr_of(&result);
                assert!(
                    stderr.contains(if fixed {
                        "Noise scale: fixed"
                    } else {
                        "no usable qualities"
                    }),
                    "{stderr}"
                );
                assert_eq!(
                    fs::read_to_string(&assignments).unwrap(),
                    expected,
                    "{count}, {fixed}, {source:?}"
                );
            }
        }
    }
}

#[test]
fn rejects_stdin_with_fastq_output_before_reading_and_requires_a_result() {
    let dir = TempDir::new().unwrap();
    let output = dir.path().join("out.fastq");
    for (source, message) in [
        ("-", "stdin requires --assignments without OUTPUT"),
        ("/dev/stdin", "regular file"),
    ] {
        Command::cargo_bin("dedup_ont_fastq")
            .unwrap()
            .arg(source)
            .arg(&output)
            .write_stdin("malformed input")
            .assert()
            .failure()
            .stderr(predicates::str::contains(message))
            .stderr(predicates::str::contains("--assignments"));
    }
    assert!(!output.exists());
    Command::cargo_bin("dedup_ont_fastq")
        .unwrap()
        .arg("-")
        .write_stdin("malformed input")
        .assert()
        .failure()
        .stderr(predicates::str::contains("--assignments"));
}

#[test]
fn rejects_assignments_that_overwrite_redirected_stdin() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.fastq");
    let contents = "@r0\nACGT\n+\nIIII\n";
    fs::write(&input, contents).unwrap();
    for source in ["-", "/dev/stdin"] {
        let output = std::process::Command::new(env!("CARGO_BIN_EXE_dedup_ont_fastq"))
            .arg(source)
            .arg("--assignments")
            .arg(&input)
            .stdin(File::open(&input).unwrap())
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("must be different files"));
        assert_eq!(fs::read_to_string(&input).unwrap(), contents);
    }
}
