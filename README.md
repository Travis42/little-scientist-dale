# DALE: Dual-seed Algorithm for Latent Enumeration

**Fast, accurate de novo DNA motif discovery for ChIP-seq data.**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Travis42/little-scientist-dale/blob/main/head_to_head_benchmark.ipynb)

**Try it in your browser** — click the badge above to open in Google Colab, then Runtime → Run All. First run takes ~3 minutes (MEME Suite install + compile). The benchmark itself runs in ~1 minute.

📄 **Paper:** [doi.org/10.5281/zenodo.21907349](https://doi.org/10.5281/zenodo.21907349)

DALE is a motif discovery algorithm written in C that outperforms STREME (the default in the MEME Suite) across 132 ENCODE transcription factors while running 11× faster. It was discovered autonomously by an LLM-based agent operating within the [Little Scientist framework](https://github.com/travis42/little-scientist-delta-v) — the agent wrote the algorithm from scratch through iterative hypothesis testing.

Distributed as a single 906 KB statically-linked binary — no dependencies, no compilation, no container required. Source code included for full reproducibility.

---

## Results

Benchmarked on **132 ENCODE K562 ChIP-seq** transcription factor datasets against STREME 5.5.5 and MEME 5.5.5.

### Accuracy

- **Shuffled negatives (132 TFs):** DALE 0.842 vs STREME 0.803 (Wilcoxon p = 4.1×10⁻⁷)
- **Genomic negatives (53 TFs):** DALE 0.891 vs STREME 0.863 (Wilcoxon p = 6.5×10⁻⁴)

### Speed

- **DALE:** 0.3s/TF (~45s for all 132 TFs) — 11× faster than STREME
- **STREME:** 3.3s/TF (~7 min for all 132)
- **MEME:** ~90s/TF (~3 hr for all 132)

Full scoreboard: [`results/SCOREBOARD.md`](results/SCOREBOARD.md)

---

## Quick Start

### Pre-compiled binary (no build required)

```bash
# Discover motifs in a ChIP-seq FASTA file
./bin/dale < input_sequences.fa
```

Works on any Linux x86_64 system with no dependencies.

### Build from source

```bash
cd src/
make
./benchmark --data ../example --tf SP1
```

Requires only a C99 compiler and libm.

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

Full walkthrough: [`docs/ALGORITHM_OVERVIEW.md`](docs/ALGORITHM_OVERVIEW.md)

---

## Origin

DALE was discovered by an autonomous LLM agent (GLM-4.7) operating within the Little Scientist framework. The agent was given a scoring function (AUROC on held-out ChIP-seq data), a smoke-test environment, and iterative feedback. It wrote the Python strategy through 118 iterations of hypothesis → implementation → empirical testing. The winning strategy was ported to C for performance. See [`sef/`](sef/) for the original strategy and agent prompt.

---

## Citation

Smith, Travis. (2026). The Little Scientist: Hypothesis-Driven Iterative Algorithm Discovery by LLM Agents. Zenodo. https://doi.org/10.5281/zenodo.21907349

## License

Apache License 2.0. See [LICENSE](LICENSE). Includes an explicit patent grant.
