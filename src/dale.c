# Copyright 2026 Travis Smith
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

/*
 * dale.c -- C port of best_so_far_strategy.py (Run #118, score 0.8503).
 *
 * Faithful port of the Pareto-front motif-discovery pipeline:
 *   1. Adaptive global background from the discovery sequences
 *   2. Multi-width scan over [6, 8, 10, 12, 14, 17]
 *   3. K-mer seed initialization:
 *        - 4 enriched seeds  (frequency, Hamming diversity >= w//3)
 *        - 2 consistent seeds (seq_count/sqrt(total+1), diversity >= w//2)
 *   4. Soft EM refinement (20 iters, pseudocount 0.5, conv 1e-4)
 *   5. Per-width z-calibration of total IC
 *   6. Discriminative score: 0.7*total_ic + 0.3*ic_density, z-bonus (0.2)
 *   7. Top-2 per width -> Pareto front (discriminative vs width efficiency)
 *   8. Quadratic padding penalty: width_eff = ic_density * (1 - pad^1.5)
 *
 * Two log bases are used deliberately:
 *   - log2 for information content (matches Python math.log2)
 *   - natural log for PWM probabilities (matches Python np.log)
 */
#define _GNU_SOURCE
#include "dale.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

/* ===================== Constants & tiny helpers ===================== */

static const char BASES[4] = {'A', 'C', 'G', 'T'};

static inline int b2i_c(char b) {
    switch (b) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

/* Encode a width-w k-mer from seq as 2 bits/base. Returns 1 on success. */
static int encode_kmer(const char *seq, int w, uint64_t *out) {
    uint64_t k = 0;
    for (int j = 0; j < w; j++) {
        int b = b2i_c(seq[j]);
        if (b < 0) return 0;
        k = (k << 2ULL) | (uint64_t)b;
    }
    *out = k;
    return 1;
}

static void decode_kmer(uint64_t key, int w, int *bases) {
    for (int j = 0; j < w; j++)
        bases[j] = (int)((key >> (2ULL * (w - 1 - j))) & 3ULL);
}

static int hamming_dist(const int *a, const int *b, int w) {
    int d = 0;
    for (int j = 0; j < w; j++) if (a[j] != b[j]) d++;
    return d;
}

static size_t next_pow2(size_t v) {
    if (v < 1) return 1;
    v--;
    v |= v >> 1; v |= v >> 2; v |= v >> 4;
    v |= v >> 8; v |= v >> 16;
    if (sizeof(size_t) > 4) v |= v >> 32;
    return v + 1;
}

static double clampd(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

/* ===================== K-mer hash table =====================
 *
 * Tracks, for each distinct k-mer key:
 *   - total_count : total occurrences across all sequences
 *   - seq_count   : number of distinct sequences containing it
 * The single-pass "last_seq" trick keeps seq_count exact: a k-mer's
 * seq_count is incremented only the first time it is seen in each sequence.
 */
#define KH_EMPTY ((uint64_t)-1)

typedef struct {
    uint64_t *keys;
    int      *total_count;
    int      *seq_count;
    int      *last_seq;
    size_t    mask;
    size_t    capacity;
} KMerTable;

static void kt_init(KMerTable *t, size_t capacity) {
    if (capacity < 16) capacity = 16;
    t->capacity = capacity;
    t->mask = capacity - 1;
    t->keys       = (uint64_t *)malloc(sizeof(uint64_t) * capacity);
    t->total_count= (int *)     malloc(sizeof(int)      * capacity);
    t->seq_count  = (int *)     malloc(sizeof(int)      * capacity);
    t->last_seq   = (int *)     malloc(sizeof(int)      * capacity);
    for (size_t i = 0; i < capacity; i++) t->keys[i] = KH_EMPTY;
}

static void kt_free(KMerTable *t) {
    free(t->keys); free(t->total_count);
    free(t->seq_count); free(t->last_seq);
    t->keys = NULL;
}

static inline void kt_add(KMerTable *t, uint64_t key, int seq_idx) {
    size_t i = (size_t)key & t->mask;
    while (t->keys[i] != KH_EMPTY && t->keys[i] != key)
        i = (i + 1) & t->mask;
    if (t->keys[i] == KH_EMPTY) {
        t->keys[i]        = key;
        t->total_count[i] = 1;
        t->seq_count[i]   = 1;
        t->last_seq[i]    = seq_idx;
    } else {
        t->total_count[i]++;
        if (t->last_seq[i] != seq_idx) {
            t->seq_count[i]++;
            t->last_seq[i] = seq_idx;
        }
    }
}

typedef struct { uint64_t key; int total_count; int seq_count; } KmerRow;

static int krow_total_cmp_desc(const void *a, const void *b) {
    const KmerRow *x = (const KmerRow *)a;
    const KmerRow *y = (const KmerRow *)b;
    if (x->total_count < y->total_count) return 1;
    if (x->total_count > y->total_count) return -1;
    /* deterministic tie-break: smaller key first */
    if (x->key < y->key) return -1;
    if (x->key > y->key) return 1;
    return 0;
}

static int krow_cons_cmp_desc(const void *a, const void *b, double ca, double cb) {
    if (ca < cb) return 1;
    if (ca > cb) return -1;
    return 0;
}

/* Wrapper struct carrying the consistency score for sorting. */
static int krow_cons_cmp_desc_g(const void *a, const void *b) {
    const KmerRow *x = (const KmerRow *)a;
    const KmerRow *y = (const KmerRow *)b;
    double cx = (double)x->seq_count / sqrt((double)x->total_count + 1.0);
    double cy = (double)y->seq_count / sqrt((double)y->total_count + 1.0);
    return krow_cons_cmp_desc(a, b, cx, cy);
}

/* Collect all non-empty entries into a freshly-malloc'd array. */
static KmerRow *kt_collect(const KMerTable *t, int *out_n) {
    KmerRow *arr = (KmerRow *)malloc(sizeof(KmerRow) * t->capacity);
    int n = 0;
    for (size_t i = 0; i < t->capacity; i++) {
        if (t->keys[i] != KH_EMPTY) {
            arr[n].key        = t->keys[i];
            arr[n].total_count= t->total_count[i];
            arr[n].seq_count  = t->seq_count[i];
            n++;
        }
    }
    *out_n = n;
    return arr;
}

/* ===================== Sequence matrix ===================== */

typedef struct {
    int  n_pos;
    int *mat;  /* n_pos * w ints, row-major; NULL if n_pos == 0 */
} SeqMatrix;

static void seq_to_matrix(const char *seq, int w, SeqMatrix *sm) {
    int len = (int)strlen(seq);
    int n_pos = len - w + 1;
    if (n_pos <= 0) {
        sm->n_pos = 0;
        sm->mat   = NULL;
        return;
    }
    int *mat = (int *)malloc(sizeof(int) * (size_t)n_pos * w);
    int valid = 0;
    for (int i = 0; i < n_pos; i++) {
        int ok = 1;
        for (int j = 0; j < w; j++) {
            int b = b2i_c(seq[i + j]);
            if (b < 0) { ok = 0; break; }
            mat[valid * w + j] = b;
        }
        if (ok) valid++;
    }
    if (valid == 0) {
        free(mat);
        sm->n_pos = 0;
        sm->mat   = NULL;
    } else {
        sm->n_pos = valid;
        sm->mat   = mat;
    }
}

/* ===================== PWM helpers ===================== */

static float *pwm_alloc(int w) {
    return (float *)malloc(sizeof(float) * (size_t)w * 4);
}

static float *pwm_copy(const float *src, int w) {
    float *dst = pwm_alloc(w);
    memcpy(dst, src, sizeof(float) * (size_t)w * 4);
    return dst;
}

static void pwm_normalize_nonneg(float *pwm, int w) {
    for (int i = 0; i < w; i++) {
        double s = 0;
        for (int b = 0; b < 4; b++) {
            if (pwm[i * 4 + b] < 0.0) pwm[i * 4 + b] = 0.0f;
            s += pwm[i * 4 + b];
        }
        if (s > 0) {
            for (int b = 0; b < 4; b++) pwm[i * 4 + b] = (float)(pwm[i * 4 + b] / s);
        } else {
            for (int b = 0; b < 4; b++) pwm[i * 4 + b] = 0.25f;
        }
    }
}

/* ===================== Seed PWM (0.25 + 2.68 for k-mer base) ===================== */

static float *build_seed_pwm(const int *kmer_bases, int w) {
    float *pwm = pwm_alloc(w);
    for (int j = 0; j < w; j++) {
        for (int b = 0; b < 4; b++) {
            double v = 0.25;
            if (b == kmer_bases[j]) v += 2.68;
            pwm[j * 4 + b] = (float)v;
        }
        double s = 0;
        for (int b = 0; b < 4; b++) s += pwm[j * 4 + b];
        for (int b = 0; b < 4; b++) pwm[j * 4 + b] = (float)(pwm[j * 4 + b] / s);
    }
    return pwm;
}

/* ===================== Soft EM (vectorized-equivalent) =====================
 *
 * 20 iterations max, pseudocount 0.5, convergence on max-abs-diff < 1e-4.
 * Background is position-independent: bg_log[b] = log(bg[b] + 1e-10).
 * Per sequence: softmax over site log-odds, accumulate expected counts.
 */
static void soft_em(const float *seed, int w, const SeqMatrix *mats, int n_seqs,
                    const double *bg, float *out_pwm)
{
    double bg_log[4];
    for (int b = 0; b < 4; b++) bg_log[b] = log(bg[b] + 1e-10);

    float *pwm = pwm_alloc(w);
    memcpy(pwm, seed, sizeof(float) * (size_t)w * 4);

    double *counts  = (double *)malloc(sizeof(double) * (size_t)w * 4);
    double *pwm_log = (double *)malloc(sizeof(double) * (size_t)w * 4);
    float  *old     = pwm_alloc(w);

    int max_npos = 0;
    for (int s = 0; s < n_seqs; s++)
        if (mats[s].n_pos > max_npos) max_npos = mats[s].n_pos;
    int bsz = max_npos > 0 ? max_npos : 1;
    double *ls   = (double *)malloc(sizeof(double) * bsz);  /* log scores  */
    double *post = (double *)malloc(sizeof(double) * bsz);  /* posteriors  */

    for (int iter = 0; iter < 20; iter++) {
        memcpy(old, pwm, sizeof(float) * (size_t)w * 4);
        for (int j = 0; j < w * 4; j++) counts[j] = 0.5;

        for (int j = 0; j < w; j++)
            for (int b = 0; b < 4; b++)
                pwm_log[j * 4 + b] = log((double)pwm[j * 4 + b] + 1e-10);

        for (int s = 0; s < n_seqs; s++) {
            const SeqMatrix *sm = &mats[s];
            if (sm->n_pos == 0) continue;

            double max_score = -INFINITY;
            for (int p = 0; p < sm->n_pos; p++) {
                const int *row = &sm->mat[(size_t)p * w];
                double sc = 0;
                for (int j = 0; j < w; j++)
                    sc += pwm_log[j * 4 + row[j]] - bg_log[row[j]];
                ls[p] = sc;
                if (sc > max_score) max_score = sc;
            }
            double sum = 0;
            for (int p = 0; p < sm->n_pos; p++) {
                post[p] = exp(ls[p] - max_score);
                sum += post[p];
            }
            if (sum > 1e-10) {
                for (int p = 0; p < sm->n_pos; p++) {
                    const int *row = &sm->mat[(size_t)p * w];
                    double pp = post[p] / sum;
                    for (int j = 0; j < w; j++)
                        counts[j * 4 + row[j]] += pp;
                }
            }
        }

        for (int j = 0; j < w; j++) {
            double cs = 0;
            for (int b = 0; b < 4; b++) cs += counts[j * 4 + b];
            for (int b = 0; b < 4; b++)
                pwm[j * 4 + b] = (float)(counts[j * 4 + b] / cs);
        }

        double maxdiff = 0;
        for (int j = 0; j < w * 4; j++) {
            double d = fabs((double)pwm[j] - (double)old[j]);
            if (d > maxdiff) maxdiff = d;
        }
        if (maxdiff < 1e-4) break;
    }

    memcpy(out_pwm, pwm, sizeof(float) * (size_t)w * 4);
    free(pwm); free(old); free(counts); free(pwm_log);
    free(ls); free(post);
}

/* ===================== IC / effective width =====================
 *
 * total_ic  = sum over positions of per-position IC (bits, log2).
 * eff_w     = number of positions with per-position IC >= 0.1.
 */
static void effective_width_and_ic(const float *pwm, int w, const double *bg,
                                   int *eff_w, double *total_ic)
{
    double tic = 0.0;
    int ew = 0;
    for (int i = 0; i < w; i++) {
        double pos_ic = 0.0;
        for (int b = 0; b < 4; b++) {
            double p = pwm[i * 4 + b];
            if (p > 1e-10 && bg[b] > 1e-10)
                pos_ic += p * log2(p / bg[b]);
        }
        if (pos_ic < 0) pos_ic = 0;
        tic += pos_ic;
        if (pos_ic >= 0.1) ew++;
    }
    *eff_w   = ew;
    *total_ic= tic;
}

static void consensus(const float *pwm, int w, char *buf) {
    for (int i = 0; i < w; i++) {
        int best = 0;
        for (int b = 1; b < 4; b++)
            if (pwm[i * 4 + b] > pwm[i * 4 + best]) best = b;
        buf[i] = BASES[best];
    }
    buf[w] = 0;
}

/* ===================== Seed selection ===================== */

#define N_ENRICHED_SEEDS    4
#define N_CONSISTENT_SEEDS  2

/* Greedily select up to `want` k-mers from `sorted` (already sorted desc by
 * their respective criterion) with pairwise Hamming distance >= min_ham.
 * Fills `out_keys` and returns the number selected. */
static int select_diverse(const KmerRow *sorted, int n, int want, int w,
                          int min_ham, uint64_t *out_keys)
{
    int selected = 0;
    int sel_bases[32][32];   /* up to 8 selections, w <= 17 */
    for (int i = 0; i < n && selected < want; i++) {
        int bases[32];
        decode_kmer(sorted[i].key, w, bases);
        int diverse = 1;
        for (int k = 0; k < selected; k++) {
            if (hamming_dist(bases, sel_bases[k], w) < min_ham) {
                diverse = 0;
                break;
            }
        }
        if (diverse) {
            out_keys[selected] = sorted[i].key;
            for (int j = 0; j < w; j++) sel_bases[selected][j] = bases[j];
            selected++;
        }
    }
    /* Fallback: if not enough diverse seeds, take the top of the sorted list. */
    for (int i = 0; selected < want && i < n; i++) {
        uint64_t k = sorted[i].key;
        int dup = 0;
        for (int kk = 0; kk < selected; kk++)
            if (out_keys[kk] == k) { dup = 1; break; }
        if (!dup) out_keys[selected++] = k;
    }
    return selected;
}

/* ===================== Candidate type ===================== */

typedef struct {
    int     w;
    float  *pwm;       /* owns the buffer (w*4) */
    int     eff_w;
    double  total_ic;
    uint64_t kmer_key;
    char    type;      /* 'e' enriched, 'c' consistent */
    /* computed during ranking */
    double  z_score;
    double  disc_score;
    double  width_eff;
    int     pareto_rank;
} Candidate;

/* Per-width candidate ordering is done with inline selection sort below
 * (stable enough for the small per-width candidate sets). */

/* ===================== Main entry point ===================== */

int discover_motifs(char **sequences, int n_sequences, int k_hint,
                    Motif **out_motifs, int *out_n_motifs)
{
    *out_motifs   = NULL;
    *out_n_motifs = 0;
    if (n_sequences <= 0 || sequences == NULL) return 0;

    /* ---- adaptive global background from the sequences ---- */
    double bg[4] = {0, 0, 0, 0};
    {
        double total = 0;
        for (int i = 0; i < n_sequences; i++) {
            for (const char *p = sequences[i]; *p; p++) {
                int b = b2i_c(*p);
                if (b >= 0) { bg[b] += 1.0; total += 1.0; }
            }
        }
        if (total > 0) for (int b = 0; b < 4; b++) bg[b] /= total;
        else           for (int b = 0; b < 4; b++) bg[b] = 0.25;
    }

    /* ---- widths ---- */
    int widths[8];
    int n_widths;
    if (k_hint > 0) {
        widths[0] = k_hint;
        n_widths  = 1;
    } else {
        const int ws[6] = {6, 8, 10, 12, 14, 17};
        memcpy(widths, ws, sizeof(ws));
        n_widths = 6;
    }

    int seq_len = (int)strlen(sequences[0]);

    /* ---- candidate accumulator ---- */
    int cap_cand = 64, n_cand = 0;
    Candidate *cands = (Candidate *)malloc(sizeof(Candidate) * cap_cand);
    memset(cands, 0, sizeof(Candidate) * cap_cand);

    /* ---- main width loop ---- */
    for (int wi = 0; wi < n_widths; wi++) {
        int w = widths[wi];
        if (w > seq_len) continue;

        /* sequence matrices for this width */
        SeqMatrix *mats = (SeqMatrix *)malloc(sizeof(SeqMatrix) * n_sequences);
        for (int i = 0; i < n_sequences; i++)
            seq_to_matrix(sequences[i], w, &mats[i]);

        /* ---- enumerate k-mers (one pass: total_count + seq_count) ---- */
        long total_kmers = 0;
        for (int s = 0; s < n_sequences; s++) {
            const char *seq = sequences[s];
            int len = (int)strlen(seq);
            for (int i = 0; i <= len - w; i++) total_kmers++;
        }
        size_t cap = next_pow2((size_t)(total_kmers > 0 ? total_kmers * 2 : 16));
        if (cap < 4096) cap = 4096;
        if (cap > (1u << 20)) cap = (1u << 20);

        KMerTable kt;
        kt_init(&kt, cap);
        for (int s = 0; s < n_sequences; s++) {
            const char *seq = sequences[s];
            int len = (int)strlen(seq);
            for (int i = 0; i <= len - w; i++) {
                uint64_t key;
                if (encode_kmer(&seq[i], w, &key)) kt_add(&kt, key, s);
            }
        }

        int n_unique = 0;
        KmerRow *rows = kt_collect(&kt, &n_unique);

        /* ---- enriched seeds: top by total_count, diversity w//3 ---- */
        KmerRow *by_total = (KmerRow *)malloc(sizeof(KmerRow) * n_unique);
        memcpy(by_total, rows, sizeof(KmerRow) * n_unique);
        qsort(by_total, n_unique, sizeof(KmerRow), krow_total_cmp_desc);

        uint64_t enr_keys[8];
        int n_enr = 0;
        if (n_unique > 0) {
            int min_ham = w / 3; if (min_ham < 1) min_ham = 1;
            n_enr = select_diverse(by_total, n_unique, N_ENRICHED_SEEDS, w,
                                   min_ham, enr_keys);
        }

        /* ---- consistent seeds: top by seq_count/sqrt(total+1), diversity w//2 ---- */
        KmerRow *by_cons = (KmerRow *)malloc(sizeof(KmerRow) * n_unique);
        memcpy(by_cons, rows, sizeof(KmerRow) * n_unique);
        qsort(by_cons, n_unique, sizeof(KmerRow), krow_cons_cmp_desc_g);

        uint64_t con_keys[8];
        int n_con = 0;
        if (n_unique > 0) {
            int min_ham = w / 2; if (min_ham < 1) min_ham = 1;
            n_con = select_diverse(by_cons, n_unique, N_CONSISTENT_SEEDS, w,
                                   min_ham, con_keys);
        }

        /* ---- build candidates from each seed ---- */
        for (int t = 0; t < n_enr + n_con; t++) {
            uint64_t key = (t < n_enr) ? enr_keys[t]
                                       : con_keys[t - n_enr];
            char ctype = (t < n_enr) ? 'e' : 'c';
            int bases[32];
            decode_kmer(key, w, bases);

            float *seed = build_seed_pwm(bases, w);
            float *pwm  = pwm_alloc(w);
            soft_em(seed, w, mats, n_sequences, bg, pwm);
            free(seed);

            int eff_w;
            double total_ic;
            effective_width_and_ic(pwm, w, bg, &eff_w, &total_ic);

            if (n_cand == cap_cand) {
                cap_cand *= 2;
                cands = (Candidate *)realloc(cands, sizeof(Candidate) * cap_cand);
                memset(&cands[n_cand], 0, sizeof(Candidate) * (cap_cand - n_cand));
            }
            cands[n_cand].w        = w;
            cands[n_cand].pwm      = pwm;
            cands[n_cand].eff_w    = eff_w;
            cands[n_cand].total_ic = total_ic;
            cands[n_cand].kmer_key = key;
            cands[n_cand].type     = ctype;
            cands[n_cand].z_score  = 0.0;
            cands[n_cand].disc_score = 0.0;
            cands[n_cand].width_eff  = 0.0;
            cands[n_cand].pareto_rank= 0;
            n_cand++;
        }

        free(by_total); free(by_cons); free(rows);
        kt_free(&kt);
        for (int i = 0; i < n_sequences; i++) free(mats[i].mat);
        free(mats);
    }

    /* ---- per-width z-score normalization of total_ic ---- */
    for (int wi = 0; wi < n_widths; wi++) {
        int w = widths[wi];
        double sum = 0;
        int cnt = 0;
        for (int i = 0; i < n_cand; i++) {
            if (cands[i].w == w) { sum += cands[i].total_ic; cnt++; }
        }
        if (cnt == 0) continue;
        double mean = sum / cnt;
        double var = 0;
        for (int i = 0; i < n_cand; i++)
            if (cands[i].w == w) {
                double d = cands[i].total_ic - mean; var += d * d;
            }
        double std = sqrt(var / cnt);
        if (std < 1e-10) std = 1.0;
        for (int i = 0; i < n_cand; i++)
            if (cands[i].w == w)
                cands[i].z_score = (cands[i].total_ic - mean) / std;
    }

    /* ---- discriminative score (per-width top-2 selection) ----
     * base = 0.7*total_ic + 0.3*ic_density
     * z_cal = clamp(z * (w/10)^0.3, -1, 1)
     * disc  = base * (1 + 0.2 * z_cal)
     */
    for (int i = 0; i < n_cand; i++) {
        double ic_density = cands[i].eff_w > 0
                              ? cands[i].total_ic / cands[i].eff_w : 0.0;
        double base = 0.7 * cands[i].total_ic + 0.3 * ic_density;
        double zcal = cands[i].z_score * pow(cands[i].w / 10.0, 0.3);
        zcal = clampd(zcal, -1.0, 1.0);
        cands[i].disc_score = base * (1.0 + 0.2 * zcal);
    }

    /* ---- select top-2 per width by discriminative score ---- */
    /* Build per-width index lists, sort, take top 2 into `winners`. */
    Candidate **winners = (Candidate **)malloc(sizeof(Candidate *) * (n_widths * 2 + 1));
    int n_win = 0;

    for (int wi = 0; wi < n_widths; wi++) {
        int w = widths[wi];
        /* collect pointers for this width */
        Candidate **pw = (Candidate **)malloc(sizeof(Candidate *) * n_cand);
        int nw = 0;
        for (int i = 0; i < n_cand; i++)
            if (cands[i].w == w) pw[nw++] = &cands[i];
        if (nw == 0) { free(pw); continue; }

        /* sort pointers by disc_score desc (selection-friendly) */
        for (int a = 0; a < nw - 1; a++) {
            int best = a;
            for (int b = a + 1; b < nw; b++)
                if (pw[b]->disc_score > pw[best]->disc_score) best = b;
            Candidate *tmp = pw[a]; pw[a] = pw[best]; pw[best] = tmp;
        }
        int top = nw < 2 ? nw : 2;
        for (int t = 0; t < top; t++) winners[n_win++] = pw[t];
        free(pw);
    }

    /* ---- width efficiency (quadratic padding penalty) ----
     * ic_density = total_ic / eff_w
     * pad_ratio  = (w - eff_w) / w
     * width_eff  = ic_density * (1 - pad_ratio^1.5)
     */
    for (int i = 0; i < n_win; i++) {
        Candidate *c = winners[i];
        double ic_density = c->eff_w > 0 ? c->total_ic / c->eff_w : 0.0;
        double pad = c->w > 0 ? (double)(c->w - c->eff_w) / c->w : 0.0;
        double penalty = 1.0 - pow(pad, 1.5);
        c->width_eff = ic_density * penalty;
    }

    /* ---- Pareto rank ----
     * j dominates i if disc_j > disc_i AND width_eff_j > width_eff_i.
     * rank_i = dominated_by + 1
     */
    for (int i = 0; i < n_win; i++) {
        int dominated = 0;
        for (int j = 0; j < n_win; j++) {
            if (i == j) continue;
            if (winners[j]->disc_score  > winners[i]->disc_score &&
                winners[j]->width_eff   > winners[i]->width_eff)
                dominated++;
        }
        winners[i]->pareto_rank = dominated + 1;
    }

    /* ---- order winners by (pareto_rank asc, disc_score desc) ----
     * Stable group sort: bucket by pareto rank (small range), keep disc order.
     */
    {
        int max_rank = 0;
        for (int i = 0; i < n_win; i++)
            if (winners[i]->pareto_rank > max_rank) max_rank = winners[i]->pareto_rank;
        Candidate **ordered = (Candidate **)malloc(sizeof(Candidate *) * (n_win > 0 ? n_win : 1));
        int n_ord = 0;
        for (int r = 1; r <= max_rank; r++) {
            /* gather all with this rank, in current order, then stable-sort
             * the gathered sublist by disc desc */
            Candidate **grp = (Candidate **)malloc(sizeof(Candidate *) * (n_win > 0 ? n_win : 1));
            int ng = 0;
            for (int i = 0; i < n_win; i++)
                if (winners[i]->pareto_rank == r) grp[ng++] = winners[i];
            for (int a = 0; a < ng - 1; a++) {
                int best = a;
                for (int b = a + 1; b < ng; b++)
                    if (grp[b]->disc_score > grp[best]->disc_score) best = b;
                Candidate *tmp = grp[a]; grp[a] = grp[best]; grp[best] = tmp;
            }
            for (int k = 0; k < ng; k++) ordered[n_ord++] = grp[k];
            free(grp);
        }
        memcpy(winners, ordered, sizeof(Candidate *) * n_ord);
        free(ordered);
    }

    /* ---- take top 5, skip eff_w < 1 ---- */
    int top_n = n_win < 5 ? n_win : 5;

    Motif *out = (Motif *)calloc(MOTIF_MAX_RETURN, sizeof(Motif));
    int n_out = 0;
    for (int i = 0; i < top_n && n_out < MOTIF_MAX_RETURN; i++) {
        Candidate *c = winners[i];
        if (c->eff_w < 1) continue;

        float *pwm = pwm_copy(c->pwm, c->w);
        pwm_normalize_nonneg(pwm, c->w);

        double pad_ratio = (double)(c->w - c->eff_w) / c->w;

        char cons[64];
        consensus(pwm, c->w, cons);

        out[n_out].pwm   = pwm;
        out[n_out].width = c->w;
        out[n_out].ic    = (float)(c->total_ic / c->w);
        snprintf(out[n_out].info, MOTIF_INFO_LEN,
                 "w=%d eff_w=%d ic=%.3f pad=%.2f disc=%.3f weff=%.3f rank=%d %s",
                 c->w, c->eff_w, c->total_ic, pad_ratio,
                 c->disc_score, c->width_eff, c->pareto_rank, cons);
        n_out++;
    }

    /* ---- free all candidate-owned memory ---- */
    for (int i = 0; i < n_cand; i++) free(cands[i].pwm);
    free(cands);
    free(winners);

    *out_motifs   = out;
    *out_n_motifs = n_out;
    return 0;
}

void motif_free_all(Motif *motifs, int n_motifs) {
    if (!motifs) return;
    for (int i = 0; i < n_motifs; i++) free(motifs[i].pwm);
    free(motifs);
}
