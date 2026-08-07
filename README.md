# DALE: Dual-seed Algorithm for Latent Enumeration

**Fast, accurate de novo DNA motif discovery for ChIP-seq data.**

DALE is a motif discovery algorithm written in C that outperforms STREME (the default in the MEME Suite) across 132 ENCODE transcription factors while running 11× faster. It was discovered autonomously by an LLM-based agent operating within the [SEF (Scientific Experiment Framework)](https://github.com/travis42/little-scientist-delta-v) — the agent wrote the algorithm from scratch through iterative hypothesis testing.

Distributed as a single 928 KB statically-linked binary — no dependencies, no compilation, no container required. Source code included for full reproducibility.

---

## Results

Benchmarked on **132 ENCODE K562 ChIP-seq** transcription factor datasets against STREME 5.5.5, MEME 5.5.5, and proto-motif-discover (earlier prototype).

### Accuracy

| Benchmark | DALE | proto | STREME | MEME |
|-----------|:----:|:-----:|:------:|:----:|
| **Shuffled neg.** (132 TFs) | **0.842** | 0.843 | 0.803 | — |
| **Genomic neg.** (53 TFs) | **0.891** | 0.874 | 0.863 | 0.881 |

DALE significantly outperforms STREME on both benchmarks (Wilcoxon p = 4.1×10⁻⁷ and 6.5×10⁻⁴) and matches/exceeds MEME on genomic negatives — while running 11–15× faster than STREME and ~300× faster than MEME.

### Speed

| Tool | Per TF | Full 132-TF | vs. STREME |
|------|:------:|:-----------:|:----------:|
| **DALE** | **0.3s** | **~45s** | **11×** |
| STREME | 3.3s | ~7 min | 1× |
| MEME | ~90s | ~3 hr | 0.003× |

### Statistical Significance (vs. STREME)

| Test | Shuffled (132) | Genomic (53) |
|------|:--------------:|:------------:|
| Wilcoxon p | 4.1×10⁻⁷ | 6.5×10⁻⁴ |

Full scoreboard: [`results/SCOREBOARD.md`](results/SCOREBOARD.md)

---

## Quick Start

### Option 1: Pre-compiled binary (no build required)

```bash
# Discover motifs in a ChIP-seq FASTA file
./bin/dale < input_sequences.fa

# With Markov-1 background scoring on negatives
./bin/dale --markov < input_sequences.fa
```

The binary is statically linked — works on any Linux x86_64 system with no dependencies.

### Option 2: Build from source

```bash
cd src/
make
./benchmark --data ../example --tf SP1
```

Requires only a C99 compiler and libm.

### Try it in your browser

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://github.com/travis42/little-scientist-dale/blob/main/demo_notebook.ipynb)

The demo notebook runs DALE live on sample ChIP-seq data and reproduces all benchmark figures.

---

## Algorithm

DALE uses a multi-resolution pipeline:

1. **Adaptive background** — nucleotide frequencies from input sequences
2. **Multi-width scan** — widths 6, 8, 10, 12, 14, 17
3. **K-mer seed initialization** — enriched seeds + consistent seeds with Hamming diversity
4. **Soft EM refinement** — 20 iterations, adaptive pseudocounts, background blending
5. **Per-width z-calibration** — normalizes IC across motif widths
6. **Discriminative scoring** — weighted IC + density with z-bonus
7. **Pareto-front ranking** — multi-objective selection (discriminative score vs. width efficiency)
8. **Quadratic padding penalty** — prevents over-extension of narrow motifs

Full algorithm walkthrough: [`docs/ALGORITHM_OVERVIEW.md`](docs/ALGORITHM_OVERVIEW.md)

---

## Origin

DALE was discovered by an autonomous LLM agent (GLM-4.7) operating within the SEF framework. The agent was given:
- A scoring function (AUROC on held-out ChIP-seq data)
- A smoke-test environment for local validation
- Iterative feedback (per-TF diagnostics)

The agent wrote the Python strategy through 118 iterations of hypothesis → implementation → empirical testing. The winning strategy was then ported to C for performance. The original Python strategy and agent prompt are included in [`sef/`](sef/) for reproducibility.

---

## Citation

> *Citation and DOI will be added upon publication.*

This repository is companion code for an upcoming paper. If you use DALE before the paper is published, please reference this repository.

## License

Apache License 2.0. See [LICENSE](LICENSE). Includes an explicit patent grant.

## Requirements

- **Binary:** Any Linux x86_64 (statically linked, no dependencies)
- **Build:** C99 compiler, libm
- **Notebook:** Python 3.8+, see `demo_notebook.ipynb`
