use flate2::read::GzDecoder;
use ont_sim::{TRUTH_HEADER, read_truth};
use std::{collections::HashMap, fs, io::Read, process::Command};

#[test]
fn cli_roundtrip_plain_gzip_and_evaluation() {
    let dir = tempfile::tempdir().unwrap();
    let truth_path = dir.path().join("truth.tsv");
    let mut expected = None;
    for name in ["reads.fastq", "reads.fastq.gz", "repeated.fastq.gz", ".gz"] {
        let path = dir.path().join(name);
        let status = Command::new(env!("CARGO_BIN_EXE_ont_simulate"))
            .args([
                "--seed",
                "17",
                "--n-molecules",
                "12",
                "--sub-rate",
                "0.01",
                "--ins-rate",
                "0.01",
                "--del-rate",
                "0.01",
                "--fastq",
            ])
            .arg(&path)
            .arg("--truth")
            .arg(&truth_path)
            .status()
            .unwrap();
        assert!(status.success());
        let mut fastq = String::new();
        if name.ends_with(".gz") {
            GzDecoder::new(fs::File::open(&path).unwrap())
                .read_to_string(&mut fastq)
                .unwrap();
        } else {
            fastq = fs::read_to_string(&path).unwrap();
        }
        if name == "repeated.fastq.gz" {
            assert_eq!(
                fs::read(&path).unwrap(),
                fs::read(dir.path().join("reads.fastq.gz")).unwrap()
            );
        }
        if let Some(ref first) = expected {
            assert_eq!(&fastq, first);
        } else {
            expected = Some(fastq.clone());
        }
        let truth = read_truth(fs::read(&truth_path).unwrap().as_slice()).unwrap();
        let lines: Vec<_> = fastq.lines().collect();
        assert_eq!(lines.len(), truth.len() * 4);
        let mut exemplars = HashMap::new();
        let mut assignments = "read_idx\texemplar_idx\n".to_owned();
        for (t, record) in truth.iter().zip(lines.chunks_exact(4)) {
            assert_eq!(record[0], format!("@{}", t.read_name));
            assert_eq!(record[1].len(), record[3].len());
            assert!(record[3].bytes().all(|q| (35..=83).contains(&q)));
            let mean =
                record[3].bytes().map(|q| (q - 33) as f64).sum::<f64>() / record[3].len() as f64;
            assert!((mean - t.mean_q).abs() < 1e-10);
            let exemplar = exemplars.entry(t.molecule_id).or_insert(t.read_idx);
            assignments.push_str(&format!("{}\t{}\n", t.read_idx, exemplar));
        }
        fs::write(dir.path().join("assign.tsv"), assignments).unwrap();
    }
    let result = Command::new(env!("CARGO_BIN_EXE_ont_evaluate"))
        .arg("--truth")
        .arg(&truth_path)
        .arg("--assignments")
        .arg(dir.path().join("assign.tsv"))
        .arg("--json")
        .arg(dir.path().join("metrics.json"))
        .output()
        .unwrap();
    assert!(result.status.success(), "{:?}", result);
    assert!(
        String::from_utf8(result.stdout)
            .unwrap()
            .contains("molecule: precision=1 recall=1")
    );
    let json = fs::read_to_string(dir.path().join("metrics.json")).unwrap();
    assert!(json.contains("\"true_molecules\":12"));
    fs::write(&truth_path, format!("{TRUTH_HEADER}\n")).unwrap();
    let invalid = Command::new(env!("CARGO_BIN_EXE_ont_evaluate"))
        .arg("--truth")
        .arg(&truth_path)
        .arg("--assignments")
        .arg(dir.path().join("assign.tsv"))
        .output()
        .unwrap();
    assert!(!invalid.status.success());
    let missing = Command::new(env!("CARGO_BIN_EXE_ont_simulate"))
        .args(["--ins-rate", "0", "--del-rate", "0", "--fastq"])
        .arg(dir.path().join("unused.fastq"))
        .arg("--truth")
        .arg(dir.path().join("unused.tsv"))
        .output()
        .unwrap();
    assert!(!missing.status.success());
    assert!(
        String::from_utf8(missing.stderr)
            .unwrap()
            .contains("--sub-rate")
    );
}
