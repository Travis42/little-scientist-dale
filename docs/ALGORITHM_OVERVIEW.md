# Run 118 Algorithm — Overview

*Analysis of `motif.c` (Run #118, score 0.8503) — written 2026-07-28*

---

## The Big Picture

The algorithm receives a set of DNA sequences (200 sequences, ~100bp each from ChIP-seq peaks). It needs to find the short recurring pattern (the "motif") that the transcription factor is binding to. It does this by: scanning for enriched k-mers as starting guesses, refining them into PWMs via EM, then ranking all candidates using a multi-criteria scoring system.

---

## Step-by-Step

### 1. Adaptive background (lines ~290-300)

Before doing anything, it counts the actual nucleotide frequencies in the input sequences. This matters because genomes aren't 25% A/C/G/T — they have biases (AT-rich vs GC-rich regions). The background is used later to compute whether a position in the PWM is actually informative or just reflects baseline composition.

### 2. Multi-width scan (widths 6, 8, 10, 12, 14, 17)

The algorithm doesn't know how wide the motif is. So it tries 6 different widths. For each width, it repeats the entire seed→EM→score pipeline independently. This is brute force but thorough — STREME does something similar.

### 3. K-mer enumeration (the hash table, lines ~90-180)

For each width `w`, it slides across every sequence and counts every w-length word (k-mer). The hash table tracks two things per unique k-mer:

- `total_count` — how many times it appears across all sequences
- `seq_count` — how many *distinct* sequences contain it

The `last_seq` trick ensures each k-mer only gets counted once per sequence, so `seq_count` reflects breadth of occurrence, not just raw abundance.

### 4. Seed selection (two strategies, 6 seeds per width)

This is the clever part. It picks 6 starting guesses per width using two different criteria:

**4a. Enriched seeds (4 per width)** — sorted by `total_count` (raw abundance), filtered so no two selected seeds are too similar (Hamming distance ≥ w/3). These capture k-mers that appear frequently — the "obvious" overrepresented patterns.

**4b. Consistent seeds (2 per width)** — sorted by `seq_count / sqrt(total_count + 1)`, filtered with stricter diversity (Hamming distance ≥ w/2). This is a different signal: it favors k-mers that appear in *many different sequences* even if they're not the most abundant overall. A k-mer that appears 50 times in one sequence but only once in 3 others gets a high total_count but low consistency score. A k-mer that appears once in 30 different sequences gets a high consistency score.

**Why both?** Enriched seeds catch strong signals. Consistent seeds catch diffuse signals that might be the real motif but drowned out by noise.

### 5. Seed → PWM (build_seed_pwm, line ~230)

Each selected k-mer becomes a starting PWM. The seed PWM gives the k-mer's bases a strong weight: starting probability is 0.25 (uniform) + 2.68 for the matching base, then normalized. So the k-mer base gets ~76% probability and each other base gets ~8%. This is a "soft" initialization — it starts strongly biased toward the seed but doesn't hard-commit.

### 6. Soft EM refinement (the core, lines ~250-310)

This is the Expectation-Maximization loop that turns the seed into a real PWM. 20 iterations max:

**E-step:** For each sequence, the algorithm computes a log-odds score for the motif starting at every possible position. High score = "the motif probably starts here." It converts these to probabilities (softmax) — essentially asking "given the current PWM, where in this sequence does the motif most likely live?"

**M-step:** It updates the PWM by accumulating weighted contributions from all positions across all sequences. Positions that scored high contribute more to the PWM. A pseudocount of 0.5 prevents any base from going to zero probability (which would make the PWM unable to recover if it made a wrong assignment early).

**Convergence:** If the PWM barely changes between iterations (max difference < 1e-4), it stops early.

The background probabilities feed into the log-odds scoring: `pwm_log[base] - bg_log[base]`. This means "how much more likely is this base under the PWM than under background?" — which is exactly what information content measures.

### 7. Information content and effective width (line ~320)

After EM converges, it computes:

- **Per-position IC:** for each position in the PWM, how many bits of information does it contribute? (KL divergence from background, in log2). A position that's 100% A when background is 25% A has IC = 2 bits. A position that matches background has IC = 0.
- **Total IC:** sum across all positions
- **Effective width (`eff_w`):** count of positions with IC ≥ 0.1 — how many positions are actually "doing work" vs just being padding

This is key: a width-17 PWM might only have 12 informative positions. The other 5 are noise. `eff_w` captures this.

### 8. Z-score normalization (per width, line ~400)

Within each width, it computes the mean and standard deviation of total IC across all 6 candidates. Each candidate gets a z-score: how many standard deviations above/below the mean is its IC?

This normalizes across widths — a width-6 PWM naturally has lower total IC than a width-17 PWM (fewer positions to accumulate information). The z-score lets the algorithm compare "this candidate is good *for its width*" rather than favoring wide PWMs.

### 9. Discriminative scoring (the ranking formula, line ~420)

Each candidate gets a composite score:

```
base = 0.7 × total_ic + 0.3 × ic_density
z_bonus = clamp(z × (w/10)^0.3, -1, +1)
disc_score = base × (1 + 0.2 × z_bonus)
```

- `ic_density` = total_ic / eff_w — information per meaningful position
- The 0.7/0.3 weighting balances raw information content against compactness
- The z-bonus (up to ±20%) rewards candidates that stand out within their width group
- The `(w/10)^0.3` term gives wider motifs slightly more z-bonus leverage

### 10. Top-2 per width → Pareto front (line ~440)

It takes the 2 best candidates from each width (so up to 12 total). Then it computes two things for each:

- `disc_score` — how discriminative is this PWM?
- `width_eff` = `ic_density × (1 - padding_ratio^1.5)` — how efficient is the width? Padding (positions that aren't informative) is penalized quadratically.

Then it does **Pareto ranking**: candidate A dominates candidate B if A is better on *both* disc_score and width_eff. Each candidate's rank = number of candidates that dominate it + 1. Rank 1 = the Pareto-optimal front.

This is a multi-objective optimization approach. Instead of collapsing everything into one score (which would lose information), it keeps two dimensions and finds the trade-off surface.

### 11. Final selection (line ~500)

Order by Pareto rank (lower = better), break ties by disc_score. Take the top 5. Skip any with eff_w < 1 (meaningless PWMs). Return them with their PWMs and metadata.

---

## What Makes This Algorithm Work

Three things stand out:

1. **Dual seed strategy.** Most motif finders use either abundance (MEME's approach) or consistency. Using both catches different kinds of motifs — strong sharp ones and diffuse broad ones.

2. **Pareto front selection.** Instead of trusting a single scoring formula to pick the best motif, it maintains multiple criteria and lets the Pareto front find the natural trade-offs. This avoids the "one metric dominates and everything else gets ignored" problem.

3. **Effective width awareness.** The algorithm doesn't just try different widths — it explicitly measures how much of each width is informative and penalizes padding. This is why it settled on w=17 for most TFs: the wide window lets EM find the real motif positions, and the padding penalty ensures the ranking doesn't reward the extra noise.

---

## Performance

The whole thing runs in ~6.5 seconds per TF in C. For comparison, STREME does ~3.3s and MEME does ~30s+ on the same data.
