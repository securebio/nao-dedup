use nao_dedup::{DedupContext, DedupParams, MinimizerParams};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashSet;

// ============================================================================
// Helper Functions
// ============================================================================

/// Generate a random DNA sequence using an existing RNG.
fn random_seq_with_rng(length: usize, rng: &mut StdRng) -> String {
    let bases = ['A', 'C', 'G', 'T'];
    (0..length)
        .map(|_| bases[rng.gen_range(0..4)])
        .collect()
}

/// Mutate a single base in a sequence at the given position.
fn mutate_base(seq: &str, pos: usize, new_base: char) -> String {
    let mut chars: Vec<char> = seq.chars().collect();
    if pos < chars.len() {
        chars[pos] = new_base;
    }
    chars.into_iter().collect()
}

/// Reverse complement a DNA sequence.
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => c,
        })
        .collect()
}

/// Helper struct to hold read data for testing
struct TestRead {
    fwd_seq: String,
    rev_seq: String,
    fwd_qual: String,
    rev_qual: String,
}

/// Process reads through DedupContext and return exemplar indices
fn run_dedup(
    reads: Vec<TestRead>,
    dedup_params: Option<DedupParams>,
    minimizer_params: Option<MinimizerParams>,
) -> Vec<usize> {
    let dedup_params = dedup_params.unwrap_or_default();
    let minimizer_params = minimizer_params.unwrap_or_default();

    let mut ctx = DedupContext::new(dedup_params, minimizer_params);

    for (idx, read) in reads.into_iter().enumerate() {
        ctx.process_read_by_index(
            idx,
            read.fwd_seq,
            read.rev_seq,
            read.fwd_qual,
            read.rev_qual,
        );
    }

    ctx.finalize();

    // Get exemplar index for each read
    (0..ctx.stats().0)
        .map(|idx| ctx.get_exemplar_index(idx))
        .collect()
}

// ============================================================================
// Basic Functionality Tests
// ============================================================================

#[test]
fn test_empty_input() {
    let reads: Vec<TestRead> = vec![];
    let exemplars = run_dedup(reads, None, None);
    assert!(exemplars.is_empty());
}

#[test]
fn test_single_read() {
    let reads = vec![TestRead {
        fwd_seq: "AAAA".to_string(),
        rev_seq: "TTTT".to_string(),
        fwd_qual: "IIII".to_string(),
        rev_qual: "IIII".to_string(),
    }];

    let exemplars = run_dedup(reads, None, None);

    assert_eq!(exemplars.len(), 1);
    assert_eq!(exemplars[0], 0); // Maps to itself
}

#[test]
fn test_identical_reads() {
    let mut rng = StdRng::seed_from_u64(42);
    let seq_f = random_seq_with_rng(100, &mut rng);
    let seq_r = random_seq_with_rng(100, &mut rng);
    let qual = "I".repeat(100);

    let reads = vec![
        TestRead {
            fwd_seq: seq_f.clone(),
            rev_seq: seq_r.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: seq_f.clone(),
            rev_seq: seq_r.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: seq_f,
            rev_seq: seq_r,
            fwd_qual: qual.clone(),
            rev_qual: qual,
        },
    ];

    let exemplars = run_dedup(reads, None, None);

    // All should map to the same exemplar
    let unique_exemplars: HashSet<_> = exemplars.iter().collect();
    assert_eq!(unique_exemplars.len(), 1);
}

#[test]
fn test_no_duplicates() {
    let mut rng = StdRng::seed_from_u64(42);
    let qual = "I".repeat(100);

    let reads = vec![
        TestRead {
            fwd_seq: random_seq_with_rng(100, &mut rng),
            rev_seq: random_seq_with_rng(100, &mut rng),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: random_seq_with_rng(100, &mut rng),
            rev_seq: random_seq_with_rng(100, &mut rng),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: random_seq_with_rng(100, &mut rng),
            rev_seq: random_seq_with_rng(100, &mut rng),
            fwd_qual: qual.clone(),
            rev_qual: qual,
        },
    ];

    let exemplars = run_dedup(reads, None, None);

    // Each should be its own exemplar
    for (idx, &exemplar_idx) in exemplars.iter().enumerate() {
        assert_eq!(idx, exemplar_idx);
    }
}

#[test]
fn test_multiple_small_clusters() {
    let mut rng = StdRng::seed_from_u64(42);

    let seq1_f = random_seq_with_rng(100, &mut rng);
    let seq1_r = random_seq_with_rng(100, &mut rng);
    let seq2_f = random_seq_with_rng(100, &mut rng);
    let seq2_r = random_seq_with_rng(100, &mut rng);
    let seq3_f = random_seq_with_rng(100, &mut rng);
    let seq3_r = random_seq_with_rng(100, &mut rng);

    let reads = vec![
        TestRead {
            fwd_seq: seq1_f.clone(),
            rev_seq: seq1_r.clone(),
            fwd_qual: "I".repeat(100), // Lower quality (Q=40)
            rev_qual: "J".repeat(100), // Lower quality (Q=41)
        },
        TestRead {
            fwd_seq: seq1_f,
            rev_seq: seq1_r,
            fwd_qual: "K".repeat(100), // Higher quality (Q=42)
            rev_qual: "K".repeat(100),
        },
        TestRead {
            fwd_seq: seq2_f,
            rev_seq: seq2_r,
            fwd_qual: "I".repeat(100), // Singleton
            rev_qual: "I".repeat(100),
        },
        TestRead {
            fwd_seq: seq3_f.clone(),
            rev_seq: seq3_r.clone(),
            fwd_qual: "I".repeat(100), // Cluster 2
            rev_qual: "I".repeat(100),
        },
        TestRead {
            fwd_seq: seq3_f,
            rev_seq: seq3_r,
            fwd_qual: "I".repeat(100), // Cluster 2
            rev_qual: "I".repeat(100),
        },
    ];

    let exemplars = run_dedup(reads, None, None);

    let unique_exemplars: HashSet<_> = exemplars.iter().collect();
    assert_eq!(unique_exemplars.len(), 3);

    // read0 and read1 should cluster together, read1 has higher quality
    assert_eq!(exemplars[0], 1); // read0 -> read1
    assert_eq!(exemplars[1], 1); // read1 -> read1

    // read2 is a singleton
    assert_eq!(exemplars[2], 2);

    // read3 and read4 should cluster together
    assert_eq!(exemplars[3], exemplars[4]);
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_realistic_reads_with_errors() {
    let mut rng = StdRng::seed_from_u64(42);

    let base_seq_f = random_seq_with_rng(150, &mut rng);
    let base_seq_r = random_seq_with_rng(150, &mut rng);
    let qual = "I".repeat(150);

    // Create reads with small differences
    let sub_f = if base_seq_f.chars().nth(50).unwrap() != 'G' {
        'G'
    } else {
        'A'
    };
    let base_seq_f_with_error = mutate_base(&base_seq_f, 50, sub_f);

    let sub_r = if base_seq_r.chars().nth(50).unwrap() != 'G' {
        'G'
    } else {
        'A'
    };
    let base_seq_r_with_error = mutate_base(&base_seq_r, 50, sub_r);

    let reads = vec![
        TestRead {
            fwd_seq: base_seq_f.clone(),
            rev_seq: base_seq_r.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: base_seq_f_with_error,
            rev_seq: base_seq_r.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: base_seq_f,
            rev_seq: base_seq_r_with_error,
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: random_seq_with_rng(150, &mut rng),
            rev_seq: random_seq_with_rng(150, &mut rng),
            fwd_qual: qual.clone(),
            rev_qual: qual,
        },
    ];

    let dedup_params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.01,
    };

    let exemplars = run_dedup(reads, Some(dedup_params), None);

    // reads 0, 1, 2 should cluster together
    assert_eq!(exemplars[0], exemplars[1]);
    assert_eq!(exemplars[1], exemplars[2]);

    // read 3 is separate
    assert_eq!(exemplars[3], 3);

    let unique_exemplars: HashSet<_> = exemplars.iter().collect();
    assert_eq!(unique_exemplars.len(), 2);
}

#[test]
fn test_approximate_match_parity() {
    // Ensure Rust handles mismatches correctly
    let seq = "A".repeat(100);
    let seq_error = format!("{}{}{}", "A".repeat(50), "T", "A".repeat(49)); // 1 mismatch (1% error)

    let reads = vec![
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "I".repeat(100),
            rev_qual: "I".repeat(100),
        },
        TestRead {
            fwd_seq: seq_error,
            rev_seq: seq,
            fwd_qual: "I".repeat(100),
            rev_qual: "I".repeat(100),
        },
    ];

    // With 2% error threshold, should match
    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.02,
    };

    let exemplars = run_dedup(reads, Some(params), None);

    // Should be clustered together
    assert_eq!(exemplars[0], exemplars[1]);
}

#[test]
fn test_approximate_match_threshold() {
    // Ensure Rust rejects matches above error threshold
    let seq = "A".repeat(100);
    let seq_error = format!("{}{}", "A".repeat(97), "TTT"); // 3 mismatches (3% error)

    let reads = vec![
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "I".repeat(100),
            rev_qual: "I".repeat(100),
        },
        TestRead {
            fwd_seq: seq_error,
            rev_seq: seq,
            fwd_qual: "I".repeat(100),
            rev_qual: "I".repeat(100),
        },
    ];

    // With 2% threshold, should NOT match
    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.02,
    };

    let exemplars = run_dedup(reads, Some(params), None);

    // Should be separate clusters
    assert_eq!(exemplars[0], 0);
    assert_eq!(exemplars[1], 1);
}

#[test]
fn test_offset_alignment_left_shift() {
    // Test that sequences match when one is shifted left by 1 base
    let seq1 = format!("G{}", "A".repeat(99));
    let seq2 = "A".repeat(99);
    let common = "T".repeat(99);

    let reads = vec![
        TestRead {
            fwd_seq: seq1,
            rev_seq: common.clone(),
            fwd_qual: "I".repeat(100),
            rev_qual: "I".repeat(99),
        },
        TestRead {
            fwd_seq: seq2,
            rev_seq: common,
            fwd_qual: "I".repeat(99),
            rev_qual: "I".repeat(99),
        },
    ];

    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.02,
    };

    let exemplars = run_dedup(reads, Some(params), None);

    // Should be clustered together
    assert_eq!(exemplars[0], exemplars[1]);
}

#[test]
fn test_offset_alignment_right_shift() {
    // Test that sequences match when one is shifted right by 1 base
    let seq1 = "A".repeat(99);
    let seq2 = format!("G{}", "A".repeat(99));
    let common = "T".repeat(99);

    let reads = vec![
        TestRead {
            fwd_seq: seq1,
            rev_seq: common.clone(),
            fwd_qual: "I".repeat(99),
            rev_qual: "I".repeat(99),
        },
        TestRead {
            fwd_seq: seq2,
            rev_seq: common,
            fwd_qual: "I".repeat(100),
            rev_qual: "I".repeat(99),
        },
    ];

    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.02,
    };

    let exemplars = run_dedup(reads, Some(params), None);

    // Should be clustered together
    assert_eq!(exemplars[0], exemplars[1]);
}

#[test]
fn test_different_length_no_match_beyond_offset() {
    // Test that sequences with length difference > max_offset don't match
    let seq1_fwd = format!("GG{}", "A".repeat(98));
    let seq2_fwd = "A".repeat(98);
    let seq1_rev = format!("TT{}", "C".repeat(98));
    let seq2_rev = "C".repeat(98);

    let reads = vec![
        TestRead {
            fwd_seq: seq1_fwd,
            rev_seq: seq1_rev,
            fwd_qual: "I".repeat(100),
            rev_qual: "I".repeat(100),
        },
        TestRead {
            fwd_seq: seq2_fwd,
            rev_seq: seq2_rev,
            fwd_qual: "I".repeat(98),
            rev_qual: "I".repeat(98),
        },
    ];

    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.01,
    };

    let exemplars = run_dedup(reads, Some(params), None);

    // Should be separate clusters
    assert_eq!(exemplars[0], 0);
    assert_eq!(exemplars[1], 1);
}

#[test]
fn test_error_fraction_with_offset() {
    // Test that offset counts toward error budget
    let seq1 = "A".repeat(200);
    let seq2 = format!("G{}", "A".repeat(199));

    // With max_offset=1 and max_error_frac=0.004:
    // Offset=-1: overlap=199, offset counts as 1, total = 1/199 ≈ 0.00503
    // Should NOT match (0.00503 > 0.004)
    let params1 = DedupParams {
        max_offset: 1,
        max_error_frac: 0.004,
    };

    let reads1 = vec![
        TestRead {
            fwd_seq: seq1.clone(),
            rev_seq: seq1.clone(),
            fwd_qual: "I".repeat(200),
            rev_qual: "I".repeat(200),
        },
        TestRead {
            fwd_seq: seq2.clone(),
            rev_seq: seq1.clone(),
            fwd_qual: "I".repeat(200),
            rev_qual: "I".repeat(200),
        },
    ];

    let exemplars1 = run_dedup(reads1, Some(params1), None);

    // Should be separate clusters
    assert_eq!(exemplars1[0], 0);
    assert_eq!(exemplars1[1], 1);

    // But with 0.006 threshold, should match (0.00503 <= 0.006)
    let params2 = DedupParams {
        max_offset: 1,
        max_error_frac: 0.006,
    };

    let reads2 = vec![
        TestRead {
            fwd_seq: seq1.clone(),
            rev_seq: seq1.clone(),
            fwd_qual: "I".repeat(200),
            rev_qual: "I".repeat(200),
        },
        TestRead {
            fwd_seq: seq2,
            rev_seq: seq1,
            fwd_qual: "I".repeat(200),
            rev_qual: "I".repeat(200),
        },
    ];

    let exemplars2 = run_dedup(reads2, Some(params2), None);

    // Should be clustered together
    assert_eq!(exemplars2[0], exemplars2[1]);
}

// ============================================================================
// Edge Cases
// ============================================================================

#[test]
fn test_windows_with_all_ns() {
    // Test that sequences with windows containing all N's are handled correctly
    let seq_with_ns = format!("{}{}{}", "A".repeat(30), "N".repeat(30), "C".repeat(90));

    let reads = vec![
        TestRead {
            fwd_seq: seq_with_ns.clone(),
            rev_seq: seq_with_ns.clone(),
            fwd_qual: "I".repeat(150),
            rev_qual: "I".repeat(150),
        },
        TestRead {
            fwd_seq: seq_with_ns.clone(),
            rev_seq: seq_with_ns,
            fwd_qual: "I".repeat(150),
            rev_qual: "I".repeat(150),
        },
    ];

    let exemplars = run_dedup(reads, None, None);

    // Should be clustered together despite N windows
    assert_eq!(exemplars[0], exemplars[1]);
}

#[test]
fn test_all_windows_ns() {
    // Test sequences where ALL windows contain only N's
    let all_ns = "N".repeat(150);
    let normal_seq = "A".repeat(150);

    let reads = vec![
        TestRead {
            fwd_seq: all_ns.clone(),
            rev_seq: all_ns.clone(),
            fwd_qual: "I".repeat(150),
            rev_qual: "I".repeat(150),
        },
        TestRead {
            fwd_seq: all_ns.clone(),
            rev_seq: all_ns,
            fwd_qual: "I".repeat(150),
            rev_qual: "I".repeat(150),
        },
        TestRead {
            fwd_seq: normal_seq.clone(),
            rev_seq: normal_seq,
            fwd_qual: "I".repeat(150),
            rev_qual: "I".repeat(150),
        },
    ];

    let exemplars = run_dedup(reads, None, None);

    // read2 (normal) should be its own cluster
    assert_eq!(exemplars[2], 2);

    // All-N reads produce no valid comparisons, so each is separate
    assert_eq!(exemplars[0], 0);
    assert_eq!(exemplars[1], 1);
}

// ============================================================================
// Advanced Scenarios
// ============================================================================

#[test]
fn test_cluster_lookup_after_leader_update() {
    // Test that cluster lookups work correctly after the best_read_id is updated
    let seq = "A".repeat(150);

    let reads = vec![
        // Read 0: lowest quality - will be the initial exemplar (Q=0)
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "!".repeat(150), // Q=0
            rev_qual: "!".repeat(150),
        },
        // Read 1: highest quality - will become the best exemplar (Q=40)
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "I".repeat(150), // Q=40
            rev_qual: "I".repeat(150),
        },
        // Read 2: medium quality - should find the existing cluster (Q=20)
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq,
            fwd_qual: "5".repeat(150), // Q=20
            rev_qual: "5".repeat(150),
        },
    ];

    let exemplars = run_dedup(reads, None, None);

    // All three should map to the same cluster (read1 has highest quality)
    let unique_exemplars: HashSet<_> = exemplars.iter().collect();
    assert_eq!(unique_exemplars.len(), 1);

    // The best exemplar should be read1 (highest quality)
    assert_eq!(exemplars[0], 1);
    assert_eq!(exemplars[1], 1);
    assert_eq!(exemplars[2], 1);
}

#[test]
fn test_cluster_lookup_multiple_updates() {
    // Test cluster lookups with multiple leader updates
    let seq = "G".repeat(150);

    let reads = vec![
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "!".repeat(150), // Q=0
            rev_qual: "!".repeat(150),
        },
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "#".repeat(150), // Q=2
            rev_qual: "#".repeat(150),
        },
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "I".repeat(150), // Q=40 - best
            rev_qual: "I".repeat(150),
        },
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "5".repeat(150), // Q=20
            rev_qual: "5".repeat(150),
        },
        TestRead {
            fwd_seq: seq.clone(),
            rev_seq: seq,
            fwd_qual: "(".repeat(150), // Q=7
            rev_qual: "(".repeat(150),
        },
    ];

    let exemplars = run_dedup(reads, None, None);

    // All should map to the same cluster
    let unique_exemplars: HashSet<_> = exemplars.iter().collect();
    assert_eq!(unique_exemplars.len(), 1);

    // All should map to read2 (highest quality)
    for (idx, &exemplar_idx) in exemplars.iter().enumerate() {
        assert_eq!(
            exemplar_idx, 2,
            "read{} mapped to {}, expected 2",
            idx, exemplar_idx
        );
    }
}

#[test]
fn test_best_read_selection_quality_length_tiebreaker() {
    // Test that when reads have the same quality, the longer one is chosen
    let qual_char = "?"; // Q=30 (Phred+33)

    // Shorter read: 100 bases
    let seq_short = "A".repeat(100);

    // Longer read: 150 bases
    let seq_long = "A".repeat(150);

    let reads = vec![
        TestRead {
            fwd_seq: seq_short.clone(),
            rev_seq: seq_short,
            fwd_qual: qual_char.repeat(100),
            rev_qual: qual_char.repeat(100),
        },
        TestRead {
            fwd_seq: seq_long.clone(),
            rev_seq: seq_long,
            fwd_qual: qual_char.repeat(150),
            rev_qual: qual_char.repeat(150),
        },
    ];

    let exemplars = run_dedup(reads, None, None);

    // Both should cluster together
    assert_eq!(exemplars[0], exemplars[1]);

    // The longer read (read1) should be chosen as the exemplar
    assert_eq!(exemplars[0], 1);
    assert_eq!(exemplars[1], 1);
}

#[test]
fn test_adapter_orientation_swapped_deduplication() {
    // Verify that the same DNA fragment with adapters in opposite orientations
    // is correctly identified as a duplicate.
    let insert_top = format!("{}{}", "GATTACA".repeat(21), "GAT"); // 150bp total
    let insert_bottom = reverse_complement(&insert_top);

    let qual = "I".repeat(150);

    // Orientation alpha
    let fwd_alpha = insert_top.clone();
    let rev_alpha = insert_bottom.clone();

    // Orientation beta (swapped)
    let fwd_beta = insert_bottom;
    let rev_beta = insert_top;

    let reads = vec![
        TestRead {
            fwd_seq: fwd_alpha.clone(),
            rev_seq: rev_alpha.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: fwd_beta.clone(),
            rev_seq: rev_beta.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual,
        },
    ];

    // Verify they're actually swapped
    assert_eq!(fwd_alpha, rev_beta);
    assert_eq!(rev_alpha, fwd_beta);

    let exemplars = run_dedup(reads, None, None);

    // They should be detected as duplicates (tolerant mode)
    assert_eq!(exemplars[0], exemplars[1]);
}

#[test]
fn test_windowing_strategy_beginning_of_read() {
    // Verify that windowing strategy anchors to the beginning of the read
    let m_params = MinimizerParams::new(7, 25, 3).unwrap();
    let d_params = DedupParams {
        max_offset: 0,
        max_error_frac: 0.07,
    };

    let read_len = 150;
    let base_seq = "C".repeat(read_len);
    let qual = "I".repeat(read_len);

    // --- Case 1: Both miss (tail-only similarity) ---
    let mut read2_tail_only: Vec<char> = base_seq.chars().collect();
    for idx in [6, 13, 20, 31, 38, 45, 56, 63, 70] {
        read2_tail_only[idx] = 'T';
    }
    let read2_tail_only: String = read2_tail_only.into_iter().collect();

    let reads1 = vec![
        TestRead {
            fwd_seq: base_seq.clone(),
            rev_seq: base_seq.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: read2_tail_only.clone(),
            rev_seq: read2_tail_only,
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
    ];

    let exemplars1 = run_dedup(reads1, Some(d_params.clone()), Some(m_params.clone()));

    assert_ne!(exemplars1[0], exemplars1[1], "Should miss tail-only similarity");

    // --- Case 2: Both hit (window 1 has similarity) ---
    let mut read2_mid_match: Vec<char> = base_seq.chars().collect();
    for idx in [6, 13, 20, 56, 63, 70] {
        read2_mid_match[idx] = 'T';
    }
    let read2_mid_match: String = read2_mid_match.into_iter().collect();

    let reads2 = vec![
        TestRead {
            fwd_seq: base_seq.clone(),
            rev_seq: base_seq,
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        TestRead {
            fwd_seq: read2_mid_match.clone(),
            rev_seq: read2_mid_match,
            fwd_qual: qual.clone(),
            rev_qual: qual,
        },
    ];

    let exemplars2 = run_dedup(reads2, Some(d_params), Some(m_params));

    assert_eq!(exemplars2[0], exemplars2[1], "Should detect similarity via window 1");
}

#[test]
fn test_rust_kmer_len_validation() {
    // Test that Rust validates kmer_len <= 32
    let result = MinimizerParams::new(33, 50, 3);
    assert!(result.is_err());
    let err_msg = result.unwrap_err();
    assert!(err_msg.contains("k-mer length must be <= 32"));
}
