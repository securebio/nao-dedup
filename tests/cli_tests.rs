use assert_cmd::Command;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use tempfile::NamedTempFile;

/// Helper function to write FASTQ records to a gzipped file.
///
/// Args:
///     path: Output file path
///     records: Slice of (header, sequence, quality) tuples
fn write_fastq_gz(path: &Path, records: &[(String, String, String)]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let encoder = GzEncoder::new(file, Compression::default());
    let mut writer = std::io::BufWriter::new(encoder);

    for (header, sequence, quality) in records {
        writeln!(writer, "@{}", header)?;
        writeln!(writer, "{}", sequence)?;
        writeln!(writer, "+")?;
        writeln!(writer, "{}", quality)?;
    }

    writer.flush()?;
    Ok(())
}

/// Helper function to read FASTQ records from a gzipped file.
///
/// Returns:
///     Vector of (header, sequence, quality) tuples (header without '@' prefix)
fn read_fastq_gz(path: &Path) -> std::io::Result<Vec<(String, String, String)>> {
    let file = File::open(path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let mut lines = reader.lines();
    let mut records = Vec::new();

    loop {
        let header = match lines.next() {
            Some(Ok(line)) => line,
            Some(Err(e)) => return Err(e),
            None => break,
        };

        let sequence = lines
            .next()
            .ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Missing sequence")
            })??;

        let _plus = lines
            .next()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Missing +"))??;

        let quality = lines
            .next()
            .ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Missing quality")
            })??;

        // Remove '@' from header
        let header_stripped = header.strip_prefix('@').unwrap_or(&header).to_string();
        records.push((header_stripped, sequence, quality));
    }

    Ok(records)
}

/// Helper function to generate random DNA sequence.
fn random_seq(length: usize, rng: &mut StdRng) -> String {
    const BASES: &[u8] = b"ACGT";
    (0..length)
        .map(|_| BASES[rng.gen_range(0..4)] as char)
        .collect()
}

#[test]
fn test_binary_basic_deduplication() {
    // Test that the binary correctly deduplicates identical read pairs.

    // Create synthetic interleaved FASTQ data
    // Format: R1, R2, R1, R2, ...
    let seq1_r1 = "ACGT".repeat(25); // 100bp
    let seq1_r2 = "TGCA".repeat(25); // 100bp
    let seq2_r1 = "GGGG".repeat(25); // Different sequence
    let seq2_r2 = "CCCC".repeat(25);
    let qual = "I".repeat(100); // High quality

    // Create input: 5 read pairs
    // Pairs 0, 1, 2 are identical (cluster 1)
    // Pairs 3, 4 are identical (cluster 2)
    let records = vec![
        // Pair 0 (cluster 1)
        ("read0_R1".to_string(), seq1_r1.clone(), qual.clone()),
        ("read0_R2".to_string(), seq1_r2.clone(), qual.clone()),
        // Pair 1 (cluster 1, duplicate)
        ("read1_R1".to_string(), seq1_r1.clone(), qual.clone()),
        ("read1_R2".to_string(), seq1_r2.clone(), qual.clone()),
        // Pair 2 (cluster 1, duplicate)
        ("read2_R1".to_string(), seq1_r1.clone(), qual.clone()),
        ("read2_R2".to_string(), seq1_r2.clone(), qual.clone()),
        // Pair 3 (cluster 2)
        ("read3_R1".to_string(), seq2_r1.clone(), qual.clone()),
        ("read3_R2".to_string(), seq2_r2.clone(), qual.clone()),
        // Pair 4 (cluster 2, duplicate)
        ("read4_R1".to_string(), seq2_r1.clone(), qual.clone()),
        ("read4_R2".to_string(), seq2_r2.clone(), qual.clone()),
    ];

    let input_file = NamedTempFile::new().unwrap();
    let output_file = NamedTempFile::new().unwrap();

    write_fastq_gz(input_file.path(), &records).unwrap();

    // Run the binary
    let mut cmd = Command::new(assert_cmd::cargo_bin!("dedup_interleaved_fastq"));
    cmd.arg(input_file.path())
        .arg(output_file.path())
        .assert()
        .success();

    // Read output
    let output_records = read_fastq_gz(output_file.path()).unwrap();

    // Should have 2 exemplar pairs (4 records total: 2 R1s + 2 R2s)
    assert_eq!(
        output_records.len(),
        4,
        "Expected 4 records (2 pairs), got {}",
        output_records.len()
    );

    // Extract pairs from output (alternating R1/R2)
    let mut output_pairs = Vec::new();
    for i in (0..output_records.len()).step_by(2) {
        let r1 = &output_records[i];
        let r2 = &output_records[i + 1];
        output_pairs.push((r1.clone(), r2.clone()));
    }

    // Should have exactly 2 unique pairs
    assert_eq!(
        output_pairs.len(),
        2,
        "Expected 2 output pairs, got {}",
        output_pairs.len()
    );

    // Verify the sequences match our two unique sequences
    let mut output_seqs = output_pairs
        .iter()
        .map(|(r1, r2)| (r1.1.clone(), r2.1.clone()))
        .collect::<Vec<_>>();
    output_seqs.sort();

    let mut expected_seqs = vec![
        (seq1_r1.clone(), seq1_r2.clone()),
        (seq2_r1.clone(), seq2_r2.clone()),
    ];
    expected_seqs.sort();

    assert_eq!(
        output_seqs, expected_seqs,
        "Output sequences don't match expected unique sequences"
    );
}

#[test]
fn test_binary_all_unique() {
    // Test that the binary handles all unique reads correctly.

    let mut rng = StdRng::seed_from_u64(12345);

    // Create 3 unique read pairs
    let mut records = Vec::new();
    for i in 0..3 {
        let r1_seq = random_seq(100, &mut rng);
        let r2_seq = random_seq(100, &mut rng);
        let qual = "I".repeat(100);
        records.push((format!("read{}_R1", i), r1_seq, qual.clone()));
        records.push((format!("read{}_R2", i), r2_seq, qual));
    }

    let input_file = NamedTempFile::new().unwrap();
    let output_file = NamedTempFile::new().unwrap();

    write_fastq_gz(input_file.path(), &records).unwrap();

    // Run the binary
    let mut cmd = Command::new(assert_cmd::cargo_bin!("dedup_interleaved_fastq"));
    cmd.arg(input_file.path())
        .arg(output_file.path())
        .assert()
        .success();

    // Read output
    let output_records = read_fastq_gz(output_file.path()).unwrap();

    // All 3 pairs should be in output (6 records)
    assert_eq!(
        output_records.len(),
        6,
        "Expected 6 records (3 unique pairs), got {}",
        output_records.len()
    );
}

#[test]
fn test_binary_quality_selection() {
    // Test that the binary selects the highest quality read as exemplar.

    let seq = "ACGT".repeat(25); // 100bp

    // Create 3 identical pairs with different quality scores
    let records = vec![
        // Pair 0: low quality (Q=10, char '+')
        ("read0_R1".to_string(), seq.clone(), "+".repeat(100)),
        ("read0_R2".to_string(), seq.clone(), "+".repeat(100)),
        // Pair 1: high quality (Q=40, char 'I')
        ("read1_R1".to_string(), seq.clone(), "I".repeat(100)),
        ("read1_R2".to_string(), seq.clone(), "I".repeat(100)),
        // Pair 2: medium quality (Q=20, char '5')
        ("read2_R1".to_string(), seq.clone(), "5".repeat(100)),
        ("read2_R2".to_string(), seq.clone(), "5".repeat(100)),
    ];

    let input_file = NamedTempFile::new().unwrap();
    let output_file = NamedTempFile::new().unwrap();

    write_fastq_gz(input_file.path(), &records).unwrap();

    // Run the binary
    let mut cmd = Command::new(assert_cmd::cargo_bin!("dedup_interleaved_fastq"));
    cmd.arg(input_file.path())
        .arg(output_file.path())
        .assert()
        .success();

    // Read output
    let output_records = read_fastq_gz(output_file.path()).unwrap();

    // Should have 1 exemplar pair (2 records)
    assert_eq!(
        output_records.len(),
        2,
        "Expected 2 records (1 pair), got {}",
        output_records.len()
    );

    // The exemplar should have high quality (all 'I')
    assert_eq!(
        output_records[0].2,
        "I".repeat(100),
        "Expected highest quality read as exemplar, got quality: {}...",
        &output_records[0].2[..10.min(output_records[0].2.len())]
    );
    assert_eq!(
        output_records[1].2,
        "I".repeat(100),
        "Expected highest quality read as exemplar, got quality: {}...",
        &output_records[1].2[..10.min(output_records[1].2.len())]
    );
}

#[test]
fn test_binary_with_custom_params() {
    // Test that the binary accepts custom deduplication parameters.

    // Create two slightly different sequences
    let seq1 = "A".repeat(100);
    let seq2 = format!("{}TTT", "A".repeat(97)); // 3 mismatches (3% error)
    let qual = "I".repeat(100);

    let records = vec![
        ("read0_R1".to_string(), seq1.clone(), qual.clone()),
        ("read0_R2".to_string(), seq1.clone(), qual.clone()),
        ("read1_R1".to_string(), seq2.clone(), qual.clone()),
        ("read1_R2".to_string(), seq1.clone(), qual.clone()), // Only R1 differs
    ];

    let input_file = NamedTempFile::new().unwrap();
    let output_file1 = NamedTempFile::new().unwrap();

    write_fastq_gz(input_file.path(), &records).unwrap();

    // Run with strict error threshold (should NOT deduplicate)
    let mut cmd = Command::new(assert_cmd::cargo_bin!("dedup_interleaved_fastq"));
    cmd.arg(input_file.path())
        .arg(output_file1.path())
        .arg("--max-error-frac")
        .arg("0.02") // 2% threshold - sequences have 3% error
        .assert()
        .success();

    let output_records1 = read_fastq_gz(output_file1.path()).unwrap();

    // Should have 2 pairs (4 records) since they don't match with strict threshold
    assert_eq!(
        output_records1.len(),
        4,
        "Expected 4 records (2 pairs), got {}",
        output_records1.len()
    );

    // Run with looser error threshold (should deduplicate)
    let output_file2 = NamedTempFile::new().unwrap();
    let mut cmd2 = Command::new(assert_cmd::cargo_bin!("dedup_interleaved_fastq"));
    cmd2.arg(input_file.path())
        .arg(output_file2.path())
        .arg("--max-error-frac")
        .arg("0.04") // 4% threshold - sequences have 3% error
        .assert()
        .success();

    let output_records2 = read_fastq_gz(output_file2.path()).unwrap();

    // Should have 1 pair (2 records) since they match with looser threshold
    assert_eq!(
        output_records2.len(),
        2,
        "Expected 2 records (1 pair), got {}",
        output_records2.len()
    );
}

#[test]
fn test_binary_empty_input() {
    // Test that the binary handles empty input files gracefully.

    let input_file = NamedTempFile::new().unwrap();
    let output_file = NamedTempFile::new().unwrap();

    // Create an empty gzipped file
    write_fastq_gz(input_file.path(), &[]).unwrap();

    // Run the binary
    let mut cmd = Command::new(assert_cmd::cargo_bin!("dedup_interleaved_fastq"));
    let assert = cmd
        .arg(input_file.path())
        .arg(output_file.path())
        .assert()
        .success();

    // Verify no NaN in output
    let output = assert.get_output();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("NaN"),
        "Binary output contains NaN: {}",
        stderr
    );

    // Verify warning message for empty input
    assert!(
        stderr.contains("Warning: No reads found in input file"),
        "Expected warning about no reads in output: {}",
        stderr
    );

    // Output should be empty
    let output_records = read_fastq_gz(output_file.path()).unwrap();
    assert_eq!(
        output_records.len(),
        0,
        "Expected empty output for empty input, got {} records",
        output_records.len()
    );
}
