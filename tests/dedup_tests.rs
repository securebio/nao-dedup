use nao_dedup::{deduplicate_read_pairs, DedupParams, MinimizerParams, ReadPair};
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

// ============================================================================
// Basic Functionality Tests
// ============================================================================

#[test]
fn test_empty_input() {
    let read_pairs: Vec<ReadPair> = vec![];
    let mapping = deduplicate_read_pairs(read_pairs, None, None);
    assert!(mapping.is_empty());
}

#[test]
fn test_single_read() {
    let rp = ReadPair {
        read_id: "read1".to_string(),
        fwd_seq: "AAAA".to_string(),
        rev_seq: "TTTT".to_string(),
        fwd_qual: "IIII".to_string(),
        rev_qual: "IIII".to_string(),
    };

    let read_pairs = vec![rp];
    let mapping = deduplicate_read_pairs(read_pairs, None, None);

    assert_eq!(mapping.len(), 1);
    assert_eq!(mapping.get("read1"), Some(&"read1".to_string()));
}

#[test]
fn test_identical_reads() {
    let mut rng = StdRng::seed_from_u64(42);
    let seq_f = random_seq_with_rng(100, &mut rng);
    let seq_r = random_seq_with_rng(100, &mut rng);
    let qual = "I".repeat(100);

    let read_pairs = vec![
        ReadPair {
            read_id: "read1".to_string(),
            fwd_seq: seq_f.clone(),
            rev_seq: seq_r.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        ReadPair {
            read_id: "read2".to_string(),
            fwd_seq: seq_f.clone(),
            rev_seq: seq_r.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        ReadPair {
            read_id: "read3".to_string(),
            fwd_seq: seq_f,
            rev_seq: seq_r,
            fwd_qual: qual.clone(),
            rev_qual: qual,
        },
    ];

    let mapping = deduplicate_read_pairs(read_pairs, None, None);

    // All should map to the same exemplar
    let exemplars: HashSet<_> = mapping.values().collect();
    assert_eq!(exemplars.len(), 1);
}

#[test]
fn test_no_duplicates() {
    let mut rng = StdRng::seed_from_u64(42);
    let qual = "I".repeat(100);

    let read_pairs = vec![
        ReadPair {
            read_id: "read1".to_string(),
            fwd_seq: random_seq_with_rng(100, &mut rng),
            rev_seq: random_seq_with_rng(100, &mut rng),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        ReadPair {
            read_id: "read2".to_string(),
            fwd_seq: random_seq_with_rng(100, &mut rng),
            rev_seq: random_seq_with_rng(100, &mut rng),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        ReadPair {
            read_id: "read3".to_string(),
            fwd_seq: random_seq_with_rng(100, &mut rng),
            rev_seq: random_seq_with_rng(100, &mut rng),
            fwd_qual: qual.clone(),
            rev_qual: qual,
        },
    ];

    let mapping = deduplicate_read_pairs(read_pairs, None, None);

    // Each should be its own exemplar
    for (read_id, exemplar_id) in &mapping {
        assert_eq!(read_id, exemplar_id);
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

    let read_pairs = vec![
        ReadPair {
            read_id: "read1".to_string(),
            fwd_seq: seq1_f.clone(),
            rev_seq: seq1_r.clone(),
            fwd_qual: "I".repeat(100),  // Lower quality (Q=40)
            rev_qual: "J".repeat(100),  // Lower quality (Q=41)
        },
        ReadPair {
            read_id: "read2".to_string(),
            fwd_seq: seq1_f,
            rev_seq: seq1_r,
            fwd_qual: "K".repeat(100),  // Higher quality (Q=42)
            rev_qual: "K".repeat(100),
        },
        ReadPair {
            read_id: "read3".to_string(),
            fwd_seq: seq2_f,
            rev_seq: seq2_r,
            fwd_qual: "I".repeat(100),  // Singleton
            rev_qual: "I".repeat(100),
        },
        ReadPair {
            read_id: "read4".to_string(),
            fwd_seq: seq3_f.clone(),
            rev_seq: seq3_r.clone(),
            fwd_qual: "I".repeat(100),  // Cluster 2
            rev_qual: "I".repeat(100),
        },
        ReadPair {
            read_id: "read5".to_string(),
            fwd_seq: seq3_f,
            rev_seq: seq3_r,
            fwd_qual: "I".repeat(100),  // Cluster 2
            rev_qual: "I".repeat(100),
        },
    ];

    let mapping = deduplicate_read_pairs(read_pairs, None, None);

    let exemplars: HashSet<_> = mapping.values().collect();
    assert_eq!(exemplars.len(), 3);

    assert_eq!(mapping.keys().len(), 5);

    // read1 and read2 should cluster together, read2 has higher quality
    assert_eq!(mapping.get("read1"), Some(&"read2".to_string()));
    assert_eq!(mapping.get("read2"), Some(&"read2".to_string()));

    // read3 is a singleton
    assert_eq!(mapping.get("read3"), Some(&"read3".to_string()));

    // read4 and read5 should cluster together
    assert_eq!(mapping.get("read4"), mapping.get("read5"));
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
    let sub_f = if base_seq_f.chars().nth(50).unwrap() != 'G' { 'G' } else { 'A' };
    let base_seq_f_with_error = mutate_base(&base_seq_f, 50, sub_f);

    let sub_r = if base_seq_r.chars().nth(50).unwrap() != 'G' { 'G' } else { 'A' };
    let base_seq_r_with_error = mutate_base(&base_seq_r, 50, sub_r);

    let read_pairs = vec![
        ReadPair {
            read_id: "read1".to_string(),
            fwd_seq: base_seq_f.clone(),
            rev_seq: base_seq_r.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        ReadPair {
            read_id: "read2".to_string(),
            fwd_seq: base_seq_f_with_error,
            rev_seq: base_seq_r.clone(),
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        ReadPair {
            read_id: "read3".to_string(),
            fwd_seq: base_seq_f,
            rev_seq: base_seq_r_with_error,
            fwd_qual: qual.clone(),
            rev_qual: qual.clone(),
        },
        ReadPair {
            read_id: "read4".to_string(),
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

    let mapping = deduplicate_read_pairs(read_pairs, Some(dedup_params), None);

    let exemplar1 = mapping.get("read1").unwrap();
    assert_eq!(mapping.get("read2"), Some(exemplar1));
    assert_eq!(mapping.get("read3"), Some(exemplar1));

    assert_eq!(mapping.get("read4"), Some(&"read4".to_string()));

    let exemplars: HashSet<_> = mapping.values().collect();
    assert_eq!(exemplars.len(), 2);
}

#[test]
fn test_approximate_match_parity() {
    // Ensure Rust handles mismatches identically to Python
    let seq = "A".repeat(100);
    let seq_error = format!("{}{}{}", "A".repeat(50), "T", "A".repeat(49)); // 1 mismatch (1% error)

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: seq.clone(),
        rev_seq: seq.clone(),
        fwd_qual: "I".repeat(100),
        rev_qual: "I".repeat(100),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: seq_error,
        rev_seq: seq,
        fwd_qual: "I".repeat(100),
        rev_qual: "I".repeat(100),
    };

    // With 2% error threshold, should match
    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.02,
    };

    let mapping = deduplicate_read_pairs(vec![rp1, rp2], Some(params), None);

    // Should be clustered together
    assert_eq!(mapping.get("r1"), mapping.get("r2"));
}

#[test]
fn test_approximate_match_threshold() {
    // Ensure Rust rejects matches above error threshold
    let seq = "A".repeat(100);
    let seq_error = format!("{}{}", "A".repeat(97), "TTT"); // 3 mismatches (3% error)

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: seq.clone(),
        rev_seq: seq.clone(),
        fwd_qual: "I".repeat(100),
        rev_qual: "I".repeat(100),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: seq_error,
        rev_seq: seq,
        fwd_qual: "I".repeat(100),
        rev_qual: "I".repeat(100),
    };

    // With 2% threshold, should NOT match
    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.02,
    };

    let mapping = deduplicate_read_pairs(vec![rp1, rp2], Some(params), None);

    // Should be separate clusters
    assert_ne!(mapping.get("r1"), mapping.get("r2"));
    assert_eq!(mapping.get("r1"), Some(&"r1".to_string()));
    assert_eq!(mapping.get("r2"), Some(&"r2".to_string()));
}

#[test]
fn test_offset_alignment_left_shift() {
    // Test that sequences match when one is shifted left by 1 base
    // Sequence 1: G + 99 A's (100 bases total)
    // Sequence 2: 99 A's (99 bases)
    // With offset=1, seq1[1:] aligns with seq2[0:], giving perfect match
    let seq1 = format!("G{}", "A".repeat(99));
    let seq2 = "A".repeat(99);
    let common = "T".repeat(99); // Same reverse read so they share buckets

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: seq1,
        rev_seq: common.clone(),
        fwd_qual: "I".repeat(100),
        rev_qual: "I".repeat(99),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: seq2,
        rev_seq: common,
        fwd_qual: "I".repeat(99),
        rev_qual: "I".repeat(99),
    };

    // With max_offset=1, these should match
    // Overlap is 99 bases, offset counts as 1 error, so 1/99 ≈ 0.0101 (1.01% error)
    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.02,
    };

    let mapping = deduplicate_read_pairs(vec![rp1, rp2], Some(params), None);

    // Should be clustered together
    assert_eq!(mapping.get("r1"), mapping.get("r2"));
}

#[test]
fn test_offset_alignment_right_shift() {
    // Test that sequences match when one is shifted right by 1 base
    // Sequence 1: 99 A's
    // Sequence 2: G + 99 A's (100 bases total)
    // With offset=-1, seq1[0:] aligns with seq2[1:], giving perfect match
    let seq1 = "A".repeat(99);
    let seq2 = format!("G{}", "A".repeat(99));
    let common = "T".repeat(99); // Same reverse read so they share buckets

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: seq1,
        rev_seq: common.clone(),
        fwd_qual: "I".repeat(99),
        rev_qual: "I".repeat(99),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: seq2,
        rev_seq: common,
        fwd_qual: "I".repeat(100),
        rev_qual: "I".repeat(99),
    };

    // With max_offset=1, these should match
    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.02,
    };

    let mapping = deduplicate_read_pairs(vec![rp1, rp2], Some(params), None);

    // Should be clustered together
    assert_eq!(mapping.get("r1"), mapping.get("r2"));
}

#[test]
fn test_different_length_no_match_beyond_offset() {
    // Test that sequences with length difference > max_offset don't match
    // Forward: Sequence 1: GG + 98 A's (100 bases) vs 98 A's
    // Reverse: Sequence 1: TT + 98 C's (100 bases) vs 98 C's
    // Both have difference of 2 bases, but max_offset=1, so should not match
    let seq1_fwd = format!("GG{}", "A".repeat(98));
    let seq2_fwd = "A".repeat(98);
    let seq1_rev = format!("TT{}", "C".repeat(98));
    let seq2_rev = "C".repeat(98);

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: seq1_fwd,
        rev_seq: seq1_rev,
        fwd_qual: "I".repeat(100),
        rev_qual: "I".repeat(100),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: seq2_fwd,
        rev_seq: seq2_rev,
        fwd_qual: "I".repeat(98),
        rev_qual: "I".repeat(98),
    };

    let params = DedupParams {
        max_offset: 1,
        max_error_frac: 0.01,
    };

    let mapping = deduplicate_read_pairs(vec![rp1, rp2], Some(params), None);

    // Should be separate clusters
    assert_eq!(mapping.get("r1"), Some(&"r1".to_string()));
    assert_eq!(mapping.get("r2"), Some(&"r2".to_string()));
}

#[test]
fn test_error_fraction_with_offset() {
    // Test that offset counts toward error budget
    // Seq1: 200 A's
    // Seq2: G + 199 A's (shifted by 1)
    let seq1 = "A".repeat(200);
    let seq2 = format!("G{}", "A".repeat(199));

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: seq1.clone(),
        rev_seq: seq1.clone(),
        fwd_qual: "I".repeat(200),
        rev_qual: "I".repeat(200),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: seq2.clone(),
        rev_seq: seq1.clone(),
        fwd_qual: "I".repeat(200),
        rev_qual: "I".repeat(200),
    };

    // With max_offset=1 and max_error_frac=0.004:
    // Offset=-1: overlap=199, offset counts as 1, total = 1/199 ≈ 0.00503
    // Should NOT match (0.00503 > 0.004)
    let params1 = DedupParams {
        max_offset: 1,
        max_error_frac: 0.004,
    };

    let mapping1 = deduplicate_read_pairs(
        vec![
            ReadPair {
                read_id: "r1".to_string(),
                fwd_seq: seq1.clone(),
                rev_seq: seq1.clone(),
                fwd_qual: "I".repeat(200),
                rev_qual: "I".repeat(200),
            },
            ReadPair {
                read_id: "r2".to_string(),
                fwd_seq: seq2.clone(),
                rev_seq: seq1.clone(),
                fwd_qual: "I".repeat(200),
                rev_qual: "I".repeat(200),
            },
        ],
        Some(params1),
        None,
    );

    // Should be separate clusters
    assert_eq!(mapping1.get("r1"), Some(&"r1".to_string()));
    assert_eq!(mapping1.get("r2"), Some(&"r2".to_string()));

    // But with 0.006 threshold, should match (0.00503 <= 0.006)
    let params2 = DedupParams {
        max_offset: 1,
        max_error_frac: 0.006,
    };

    let mapping2 = deduplicate_read_pairs(vec![rp1, rp2], Some(params2), None);

    // Should be clustered together
    assert_eq!(mapping2.get("r1"), mapping2.get("r2"));
}

// ============================================================================
// Edge Cases
// ============================================================================

#[test]
fn test_windows_with_all_ns() {
    // Test that sequences with windows containing all N's are handled correctly
    // With default params: num_windows=4, window_len=25
    // Create sequences where middle window is all N's
    let seq_with_ns = format!("{}{}{}", "A".repeat(30), "N".repeat(30), "C".repeat(90));

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: seq_with_ns.clone(),
        rev_seq: seq_with_ns.clone(),
        fwd_qual: "I".repeat(150),
        rev_qual: "I".repeat(150),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: seq_with_ns.clone(),
        rev_seq: seq_with_ns,
        fwd_qual: "I".repeat(150),
        rev_qual: "I".repeat(150),
    };

    let mapping = deduplicate_read_pairs(vec![rp1, rp2], None, None);

    // Should be clustered together despite N windows
    assert_eq!(
        mapping.get("r1"),
        mapping.get("r2"),
        "Sequences with N windows should be clustered"
    );
}

#[test]
fn test_all_windows_ns() {
    // Test sequences where ALL windows contain only N's
    // Rust: All-N reads produce no valid comparisons, so each is separate
    let all_ns = "N".repeat(150);

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: all_ns.clone(),
        rev_seq: all_ns.clone(),
        fwd_qual: "I".repeat(150),
        rev_qual: "I".repeat(150),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: all_ns.clone(),
        rev_seq: all_ns.clone(),
        fwd_qual: "I".repeat(150),
        rev_qual: "I".repeat(150),
    };

    let normal_seq = "A".repeat(150);
    let rp3 = ReadPair {
        read_id: "r3".to_string(),
        fwd_seq: normal_seq.clone(),
        rev_seq: normal_seq,
        fwd_qual: "I".repeat(150),
        rev_qual: "I".repeat(150),
    };

    let mapping = deduplicate_read_pairs(vec![rp1, rp2, rp3], None, None);

    // r3 should always be its own cluster
    assert_eq!(mapping.get("r3"), Some(&"r3".to_string()));

    // Rust: all-N reads produce no valid comparisons for all-N reads
    assert_eq!(mapping.get("r1"), Some(&"r1".to_string()), "r1 should be its own exemplar");
    assert_eq!(mapping.get("r2"), Some(&"r2".to_string()), "r2 should be its own exemplar");
}

// ============================================================================
// Advanced Scenarios
// ============================================================================

#[test]
fn test_cluster_lookup_after_leader_update() {
    // Test that cluster lookups work correctly after the best_read_id is updated
    // This is a regression test for a bug where the cluster hash table uses the
    // initial exemplar ID as its key, but lookups incorrectly compared against
    // best_read_id, which can change during processing.
    //
    // Scenario:
    // 1. Read A (low quality) creates cluster keyed by "readA"
    // 2. Read B (high quality) matches A, updates best_read_id to "readB"
    // 3. Read C (low quality) matches A, should find the same cluster

    let seq = "A".repeat(150);

    // Read A: lowest quality - will be the initial exemplar (Q=0)
    let rp_a = ReadPair {
        read_id: "readA".to_string(),
        fwd_seq: seq.clone(),
        rev_seq: seq.clone(),
        fwd_qual: "!".repeat(150), // Q=0
        rev_qual: "!".repeat(150),
    };

    // Read B: highest quality - will become the best exemplar (Q=40)
    let rp_b = ReadPair {
        read_id: "readB".to_string(),
        fwd_seq: seq.clone(),
        rev_seq: seq.clone(),
        fwd_qual: "I".repeat(150), // Q=40
        rev_qual: "I".repeat(150),
    };

    // Read C: medium quality - should find the existing cluster (Q=20)
    let rp_c = ReadPair {
        read_id: "readC".to_string(),
        fwd_seq: seq.clone(),
        rev_seq: seq,
        fwd_qual: "5".repeat(150), // Q=20
        rev_qual: "5".repeat(150),
    };

    // Process in order A, B, C
    let read_pairs = vec![rp_a, rp_b, rp_c];
    let mapping = deduplicate_read_pairs(read_pairs, None, None);

    // All three should map to the same cluster
    let exemplars: HashSet<_> = mapping.values().collect();
    assert_eq!(
        exemplars.len(),
        1,
        "Expected 1 cluster, got {}: {:?}",
        exemplars.len(),
        exemplars
    );

    // The best exemplar should be readB (highest quality)
    assert_eq!(mapping.get("readA"), Some(&"readB".to_string()));
    assert_eq!(mapping.get("readB"), Some(&"readB".to_string()));
    assert_eq!(mapping.get("readC"), Some(&"readB".to_string()));
}

#[test]
fn test_cluster_lookup_multiple_updates() {
    // Test cluster lookups with multiple leader updates
    // This extends the previous test with more reads to ensure the bug
    // doesn't create multiple duplicate clusters

    let seq = "G".repeat(150);

    let read_pairs = vec![
        ReadPair {
            read_id: "read1".to_string(),
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "!".repeat(150), // Q=0
            rev_qual: "!".repeat(150),
        },
        ReadPair {
            read_id: "read2".to_string(),
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "#".repeat(150), // Q=2
            rev_qual: "#".repeat(150),
        },
        ReadPair {
            read_id: "read3".to_string(),
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "I".repeat(150), // Q=40 - best
            rev_qual: "I".repeat(150),
        },
        ReadPair {
            read_id: "read4".to_string(),
            fwd_seq: seq.clone(),
            rev_seq: seq.clone(),
            fwd_qual: "5".repeat(150), // Q=20
            rev_qual: "5".repeat(150),
        },
        ReadPair {
            read_id: "read5".to_string(),
            fwd_seq: seq.clone(),
            rev_seq: seq,
            fwd_qual: "(".repeat(150), // Q=7
            rev_qual: "(".repeat(150),
        },
    ];

    let mapping = deduplicate_read_pairs(read_pairs, None, None);

    // All should map to the same cluster
    let exemplars: HashSet<_> = mapping.values().collect();
    assert_eq!(
        exemplars.len(),
        1,
        "Expected 1 cluster, got {}: {:?}",
        exemplars.len(),
        exemplars
    );

    // All should map to read3 (highest quality)
    for read_id in &["read1", "read2", "read3", "read4", "read5"] {
        assert_eq!(
            mapping.get(*read_id),
            Some(&"read3".to_string()),
            "{} mapped to {:?}, expected read3",
            read_id,
            mapping.get(*read_id)
        );
    }
}

#[test]
fn test_best_read_selection_quality_length_tiebreaker() {
    // Test that when reads have the same quality, the longer one is chosen as exemplar
    // This test verifies the scoring formula: score = mean_quality * 1000 + length
    //
    // Scenario:
    // 1. Read A (shorter, Q=30) creates cluster
    // 2. Read B (longer, Q=30) matches A and should become the exemplar

    let qual_char = "?"; // Q=30 (Phred+33)

    // Shorter read: 100 bases
    let seq_short = "A".repeat(100);
    let rp_a = ReadPair {
        read_id: "readA".to_string(),
        fwd_seq: seq_short.clone(),
        rev_seq: seq_short,
        fwd_qual: qual_char.repeat(100),
        rev_qual: qual_char.repeat(100),
    };

    // Longer read: 150 bases
    let seq_long = "A".repeat(150);
    let rp_b = ReadPair {
        read_id: "readB".to_string(),
        fwd_seq: seq_long.clone(),
        rev_seq: seq_long,
        fwd_qual: qual_char.repeat(150),
        rev_qual: qual_char.repeat(150),
    };

    // Process shorter first, then longer
    let read_pairs = vec![rp_a, rp_b];
    let mapping = deduplicate_read_pairs(read_pairs, None, None);

    // Both should cluster together
    assert_eq!(
        mapping.get("readA"),
        mapping.get("readB"),
        "Reads should cluster together: A->{:?}, B->{:?}",
        mapping.get("readA"),
        mapping.get("readB")
    );

    // The longer read (readB) should be chosen as the exemplar
    assert_eq!(
        mapping.get("readA"),
        Some(&"readB".to_string()),
        "Expected readB (longer) as exemplar, but got {:?}",
        mapping.get("readA")
    );
    assert_eq!(mapping.get("readB"), Some(&"readB".to_string()));
}

#[test]
fn test_adapter_orientation_swapped_deduplication() {
    // Verify that the same DNA fragment with adapters in opposite orientations
    // is correctly identified as a duplicate.
    //
    // When adapters attach to a double-stranded insert in opposite orientations,
    // we get reads that are swapped but NOT reverse complemented:
    //
    // Orientation alpha (P5 on top strand):
    //     Forward read: beginning of top strand
    //     Reverse read: beginning of bottom strand (reported as-is by sequencer)
    //
    // Orientation beta (P5 on bottom strand):
    //     Forward read: beginning of bottom strand
    //     Reverse read: beginning of top strand (reported as-is by sequencer)
    //
    // This means: (F_alpha, R_alpha) should match (R_beta, F_beta)
    // with NO reverse complement needed.

    // Create a 150bp "insert" - this represents the actual DNA fragment
    // Use a pattern that's clearly directional (not palindromic)
    let insert_top = format!("{}{}", "GATTACA".repeat(21), "GAT"); // 150bp total
    let insert_bottom = reverse_complement(&insert_top);

    let qual = "I".repeat(150);

    // Orientation alpha: P5 attached to top strand
    // Forward reads from top strand, reverse reads from bottom strand
    let fwd_alpha = insert_top.clone();
    let rev_alpha = insert_bottom.clone();

    // Orientation beta: P5 attached to bottom strand (insert rotated 180°)
    // Forward reads from bottom strand, reverse reads from top strand
    let fwd_beta = insert_bottom.clone();
    let rev_beta = insert_top.clone();

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: fwd_alpha.clone(),
        rev_seq: rev_alpha.clone(),
        fwd_qual: qual.clone(),
        rev_qual: qual.clone(),
    };

    let rp2 = ReadPair {
        read_id: "r2".to_string(),
        fwd_seq: fwd_beta.clone(),
        rev_seq: rev_beta.clone(),
        fwd_qual: qual.clone(),
        rev_qual: qual,
    };

    // Verify they're actually swapped (not identical)
    assert_eq!(fwd_alpha, rev_beta);
    assert_eq!(rev_alpha, fwd_beta);
    assert_ne!(fwd_alpha, fwd_beta); // Different orientations

    let mapping = deduplicate_read_pairs(vec![rp1, rp2], None, None);

    // They should be detected as duplicates (tolerant mode)
    assert_eq!(
        mapping.get("r1"),
        mapping.get("r2"),
        "Adapter-swapped orientations not deduplicated. R1: {:?}, R2: {:?}",
        mapping.get("r1"),
        mapping.get("r2")
    );
}

#[test]
fn test_windowing_strategy_beginning_of_read() {
    // Verify that windowing strategy anchors to the beginning of the read
    //
    // Both implementations use adjacent windows starting from position 0,
    // focusing on the most stable region of the read (since quality drops off
    // as you get farther into a read, so more trimming is likely).
    //
    // With num_windows=3, window_len=25 on 150bp reads:
    // Windows: [0-25], [25-50], [50-75]
    //
    // Case 1: Similarity only in the tail [75-150] → miss
    // Case 2: Similarity in window 1 [25-50] → detect

    // Use parameters that allow 9 mismatches in 150bp (6% error)
    let m_params = MinimizerParams::new(7, 25, 3).unwrap();
    let d_params = DedupParams {
        max_offset: 0,
        max_error_frac: 0.07,
    };

    let read_len = 150;
    let base_seq = "C".repeat(read_len);
    let qual = "I".repeat(read_len);

    // --- Case 1: Both miss (tail-only similarity) ---
    // Break all k-mers in the first three windows [0-75]
    // but leave the tail [75-150] clean.
    let mut read2_tail_only: Vec<char> = base_seq.chars().collect();
    // 3 changes per 25bp window to disrupt all 7-mers
    for idx in [6, 13, 20, 31, 38, 45, 56, 63, 70] {
        read2_tail_only[idx] = 'T';
    }
    let read2_tail_only: String = read2_tail_only.into_iter().collect();

    let rp1 = ReadPair {
        read_id: "r1".to_string(),
        fwd_seq: base_seq.clone(),
        rev_seq: base_seq.clone(),
        fwd_qual: qual.clone(),
        rev_qual: qual.clone(),
    };

    let rp2_tail = ReadPair {
        read_id: "r2_tail".to_string(),
        fwd_seq: read2_tail_only.clone(),
        rev_seq: read2_tail_only,
        fwd_qual: qual.clone(),
        rev_qual: qual.clone(),
    };

    // Both implementations: All windows contain mismatches -> No shared minimizers
    let mapping = deduplicate_read_pairs(
        vec![rp1.clone(), rp2_tail],
        Some(d_params.clone()),
        Some(m_params.clone()),
    );

    assert_ne!(
        mapping.get("r1"),
        mapping.get("r2_tail"),
        "Should miss tail-only similarity (windows only cover first 75bp)"
    );

    // --- Case 2: Both hit (window 1 has similarity) ---
    // Break k-mers in windows 0 and 2, but leave window 1 [25-50] clean
    let mut read2_mid_match: Vec<char> = base_seq.chars().collect();
    for idx in [6, 13, 20, 56, 63, 70] {
        read2_mid_match[idx] = 'T';
    }
    let read2_mid_match: String = read2_mid_match.into_iter().collect();

    let rp2_mid = ReadPair {
        read_id: "r2_mid".to_string(),
        fwd_seq: read2_mid_match.clone(),
        rev_seq: read2_mid_match,
        fwd_qual: qual.clone(),
        rev_qual: qual,
    };

    // Both implementations: Window 1 [25-50] is intact -> Shares minimizer
    let mapping2 = deduplicate_read_pairs(
        vec![rp1, rp2_mid],
        Some(d_params),
        Some(m_params),
    );

    assert_eq!(
        mapping2.get("r1"),
        mapping2.get("r2_mid"),
        "Should detect similarity via shared minimizer in window 1"
    );
}

#[test]
fn test_rust_kmer_len_validation() {
    // Test that Rust validates kmer_len <= 32

    // kmer_len > 32 should return an error
    let result = MinimizerParams::new(33, 50, 3);
    assert!(result.is_err());
    let err_msg = result.unwrap_err();
    assert!(err_msg.contains("k-mer length must be <= 32"));
}
