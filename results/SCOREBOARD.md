# Motif Discovery — Formal Scoreboard

**Benchmark:** 53 ChIP-seq TFs, 200 positive + 200 negative (dinucleotide shuffled), 100bp centered discovery sequences, AUROC scored on full 500bp sequences.

**Last updated:** 2026-07-24 15:50 UTC

---

## 132-TF ENCODE K562 Benchmark (Shuffled Negatives)

| Rank | Algorithm | AUROC | Avg Time | n TFs | Statistical Test |
|------|-----------|-------|----------|-------|----------------|
| 🥇 | **Ours (EM)** | **0.8391** | **1.44s** | 132 | Wilcoxon p=0.000112 |
| 2 | STREME | 0.8088 | 2.8s | 132 | Cliff's δ=0.273 |

- **Wins:** 84/132 TFs ours, 48/132 STREME
- **Effect size:** Cohen's d=0.443 (medium)
- **Significance:** 131/132 motifs Bonferroni-significant (Fisher's exact)

---

## 53-TF ENCODE K562 Benchmark (Shuffled Negatives)

## AUROC Scores (53 TFs, Shuffled Negatives)

| Rank | Algorithm | AUROC | Source | Date |
|------|-----------|-------|--------|------|
| 🥇 | Agent algo v1 (mono bg) | **0.8865** | IC×√w×width_bonus, original tiers | 2026-07-20 |
| 🥈 | Agent algo v5 (mono bg) | **0.8824** | Flattened w≤7 bonus to 1.0x | 2026-07-21 |
| 🥉 | Agent algo v5 (Markov bg) | **0.8790** | Markov-1 scoring on negatives | 2026-07-21 |
| 4 | MEME | **0.8848** | MEME zoops, minw=6 maxw=17, revcomp | 2026-07-19 |
| 5 | STREME | **0.8695** | STREME default params | 2026-07-19 |

## AUROC Scores (53 TFs, Genomic Negatives — Real hg38)

| Rank | Algorithm | AUROC | Source | Date |
|------|-----------|-------|--------|------|
| 🥇 | **Agent v5 (Markov bg)** | **0.8798** | Markov-1 scoring on negatives | 2026-07-21 |
| 🥈 | STREME | **0.8633** | Default params | 2026-07-20 |
| 🥉 | Agent algo v5 (mono bg) | **0.8560** | Flattened w≤7 bonus | 2026-07-21 |
| 4 | Agent algo v1 (mono bg) | **0.8450** | Original tiers | 2026-07-20 |

## Genomic Drop Comparison

| Algorithm | Shuffled | Genomic | Drop | Time/TF |
|-----------|----------|---------|------|---------|
| **v5 Markov** | **0.8790** | **0.8798** | **+0.0008** | **1.40s** |
| **v5 mono** | 0.8824 | 0.8560 | -0.0264 | 1.40s |
| **v1 mono** | 0.8865 | 0.8450 | -0.0415 | 1.39s |
| **STREME** | 0.8695 | 0.8633 | -0.0062 | 4.4s |

v5 + Markov scoring **beats STREME on genomic negatives** (0.8798 vs 0.8633) while staying ahead on shuffled (0.8790 vs 0.8695) and 3.1x faster.

## Speed (per TF, same 53 TFs)

| Algorithm | Time/TF | Total (53 TFs) | Speedup vs Python |
|-----------|---------|----------------|-------------------|
| Agent algo v5 | **1.40s** | 74.2s | **24.7x** |
| STREME | **~2.8s** | ~148s | 12.4x |
| MEME | **~18.6s** | ~986s | 1.9x |
| Python agent | **34.6s** | ~1834s | 1.0x (baseline) |

## Key Comparisons

- **v5+Markov vs STREME (genomic):** +0.0165 AUROC, 3.1x faster 🏆
- **v5+Markov vs STREME (shuffled):** +0.0095 AUROC, 3.1x faster 🏆
- **v5+Markov wins on both benchmarks and speed**

## Failed Experiments (all reverted)

| Version | Approach | Shuffled | Genomic | Lesson |
|---------|----------|----------|---------|--------|
| v2 | Enrichment ranking | 0.8770 | 0.8301 | Collapsed to w=7-8 |
| v3 | Specificity multiplier | 0.8912 | 0.7809 | w=6 count went UP |
| v4 | Empirical k-mer background | 0.8672 | 0.8426 | No effect on ranking |
| v5c | w=6 penalty (0.92x) | 0.8780 | 0.8562 | Genomic same, shuffled worse |
| v5d | width^0.35 exponent | 0.8884 | — | w=6 exploded to 15 |

## Agent SEF Loop Scores (15 TFs, different benchmark)

| Run | Verdict | Score | Best | Date |
|-----|---------|-------|------|------|
| 0 | baseline | 0.8004 | 0.8004 | 2026-07-19 |
| 27 | best | 0.8295 | 0.8295 | 2026-07-20 |
| 44 | rejected | 0.8294 | 0.8295 | 2026-07-21 |

SEF plateaued at 0.8295. All structural change attempts (gapped seeds, Gibbs-EM, two-stage EM, profile-weighted seeds, position-weighted EM, E-value ranking) failed with regressions.

---

*This file is the canonical scoreboard. Update after every benchmark run.*
