# DALE — Formal Scoreboard

**Benchmark:** ENCODE K562 ChIP-seq transcription factor datasets.
AUROC scored against dinucleotide-shuffled negatives (132 TFs) and real hg38 genomic negatives (53 TFs).

---

## 132-TF Shuffled Negatives

| Rank | Algorithm | AUROC | Avg Time | n TFs | Statistical Test |
|------|-----------|-------|----------|-------|------------------|
| 🥇 | **DALE** | **0.842** | **0.3s** | 132 | Wilcoxon p=4.1×10⁻⁷ vs STREME |
| 2 | proto (prototype) | 0.843 | 1.6s | 132 | — |
| 3 | STREME 5.5.5 | 0.803 | 3.3s | 132 | — |

- **Wins:** 85/132 DALE, 43/132 STREME, 4 ties
- **Speed:** 11× faster than STREME

## 53-TF Genomic Negatives (Real hg38)

| Rank | Algorithm | AUROC | Avg Time | n TFs | Statistical Test |
|------|-----------|-------|----------|-------|------------------|
| 🥇 | **DALE** | **0.891** | **0.3s** | 53 | Wilcoxon p=6.5×10⁻⁴ vs STREME |
| 2 | MEME 5.5.5 | 0.881 | ~90s | 53 | — |
| 3 | STREME 5.5.5 | 0.863 | 3.3s | 53 | — |
| 4 | proto (prototype) | 0.874 | — | 53 | — |

- **Wins:** 35/53 DALE, 16/53 STREME
- **Speed:** 15× faster than STREME, ~300× faster than MEME

---

*Data source: `benchmark_data/stats.json`*
