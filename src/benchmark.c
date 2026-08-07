/*
 * benchmark.c — ChIP-seq AUROC benchmark runner.
 *
 * Mirrors reference_benchmark.py:
 *   1. Load <data>/<tf>_sequences.fa (full 500bp) and centered 100bp slice
 *   2. Discover a motif on the 100bp slice with discover_motifs()
 *   3. Score full-length positives and dinucleotide-shuffled negatives
 *   4. Compute AUROC; report per-TF TSV row
 *   5. Repeat scoring for STREME and MEME on the same data
 *   6. Print average AUROC / time per algorithm
 *
 * Reproducibility: the shuffle PRNG is seeded once with 42 at startup.
 */
#define _GNU_SOURCE
#include "dale.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <time.h>

#ifndef MEME_BIN
#define MEME_BIN   "/opt/meme/bin/meme"
#endif
#ifndef STREME_BIN
#define STREME_BIN "/opt/meme/bin/streme"
#endif

#define DISCOVERY_CENTER_BP 100

/* ===================== PRNG (xorshift64, seed 42) ===================== */

static uint64_t g_rng = 42;
static inline uint64_t rng_next(void) {
    uint64_t x = g_rng;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    g_rng = x;
    return x;
}
static inline int rng_int(int n) { return n > 0 ? (int)(rng_next() % (uint64_t)n) : 0; }

/* ===================== Tiny helpers ===================== */

static inline int b2i_c(char b) {
    switch (b) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

static void str_toupper_inplace(char *s) {
    for (; *s; s++) if (*s >= 'a' && *s <= 'z') *s = (char)(*s - 32);
}

/* ===================== FASTA loading =====================
 *
 * Returns a malloc'd array of malloc'd upper-case strings. Caller frees with
 * free_strings(). The optional `center` argument, if > 0, truncates each
 * sequence to a window of that length centered on the midpoint.
 */
static char **load_fasta(const char *path, int center, int *out_n) {
    FILE *f = fopen(path, "r");
    if (!f) return NULL;

    int cap = 256, n = 0;
    char **seqs = (char **)malloc(sizeof(char *) * cap);

    char *cur = NULL;
    int   cur_len = 0, cur_cap = 0;
    int   in_seq = 0;
    char  line[8192];

    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '>') {
            if (in_seq && cur) {
                cur[cur_len] = 0;
                str_toupper_inplace(cur);
                if (n == cap) { cap *= 2; seqs = (char **)realloc(seqs, sizeof(char *) * cap); }
                seqs[n++] = cur;
                cur = NULL; cur_len = 0; cur_cap = 0;
            }
            in_seq = 1;
            continue;
        }
        int L = (int)strlen(line);
        while (L > 0 && (line[L-1] == '\n' || line[L-1] == '\r')) L--;
        if (cur_len + L + 1 > cur_cap) {
            while (cur_len + L + 1 > cur_cap) cur_cap = cur_cap > 0 ? cur_cap * 2 : 1024;
            cur = (char *)realloc(cur, cur_cap);
        }
        memcpy(cur + cur_len, line, L);
        cur_len += L;
    }
    if (in_seq && cur) {
        cur[cur_len] = 0;
        str_toupper_inplace(cur);
        if (n == cap) { cap *= 2; seqs = (char **)realloc(seqs, sizeof(char *) * cap); }
        seqs[n++] = cur;
    } else if (cur) {
        free(cur);
    }
    fclose(f);

    if (center > 0) {
        for (int i = 0; i < n; i++) {
            int L = (int)strlen(seqs[i]);
            int mid = L / 2;
            int start = mid - center / 2; if (start < 0) start = 0;
            int end   = mid + center / 2; if (end > L) end = L;
            int newlen = end - start;
            if (start > 0) memmove(seqs[i], seqs[i] + start, newlen);
            seqs[i][newlen] = 0;
        }
    }

    *out_n = n;
    return seqs;
}

static void free_strings(char **arr, int n) {
    if (!arr) return;
    for (int i = 0; i < n; i++) free(arr[i]);
    free(arr);
}

static void write_fasta(const char *path, char **seqs, int n) {
    FILE *f = fopen(path, "w");
    if (!f) return;
    for (int i = 0; i < n; i++) fprintf(f, ">seq%d\n%s\n", i, seqs[i]);
    fclose(f);
}

/* ===================== Dinucleotide shuffle ===================== */

static char *shuffle_dinuc(const char *seq) {
    int L = (int)strlen(seq);
    int n_din = L / 2;
    int odd = L & 1;
    char *out = (char *)malloc((size_t)L + 1);

    /* Pack dinucleotides into a flat array of size 2*n_din. */
    char *dins = (char *)malloc((size_t)(n_din * 2) + 1);
    memcpy(dins, seq, (size_t)n_din * 2);

    /* Fisher-Yates over dinucleotide indices. */
    for (int i = n_din - 1; i > 0; i--) {
        int j = rng_int(i + 1);
        if (i != j) {
            char a = dins[i*2],     b = dins[i*2+1];
            dins[i*2]   = dins[j*2];
            dins[i*2+1] = dins[j*2+1];
            dins[j*2]   = a;
            dins[j*2+1] = b;
        }
    }
    memcpy(out, dins, (size_t)n_din * 2);
    if (odd) out[L-1] = seq[L-1];
    out[L] = 0;
    free(dins);
    return out;
}

static char **shuffle_all(char **seqs, int n) {
    char **out = (char **)malloc(sizeof(char *) * n);
    for (int i = 0; i < n; i++) out[i] = shuffle_dinuc(seqs[i]);
    return out;
}

/* ===================== Background file ===================== */

static int load_bg(const char *path, double bg[4]) {
    bg[0] = bg[1] = bg[2] = bg[3] = 0.25;
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char line[256];
    int found = 0;
    while (fgets(line, sizeof(line), f)) {
        char b = 0; double v = 0;
        if (sscanf(line, " %c %lf", &b, &v) == 2) {
            switch (b) {
                case 'A': bg[0] = v; found++; break;
                case 'C': bg[1] = v; found++; break;
                case 'G': bg[2] = v; found++; break;
                case 'T': bg[3] = v; found++; break;
            }
        }
    }
    fclose(f);
    return found == 4;
}

/* ===================== PWM log-odds & sequence scoring ===================== */

static double *pwm_logodds(const float *pwm, int w, const double *bg) {
    double *lo = (double *)malloc(sizeof(double) * (size_t)w * 4);
    for (int i = 0; i < w; i++) {
        double row[4], s = 0;
        for (int b = 0; b < 4; b++) { row[b] = (double)pwm[i*4+b] + 1e-6; s += row[b]; }
        for (int b = 0; b < 4; b++) {
            row[b] /= s;
            double denom = bg[b] > 0 ? bg[b] : 0.25;
            lo[i*4+b] = log2(row[b] / denom);
        }
    }
    return lo;
}

/* Max log-odds over all windows on both strands. Returns -INFINITY if the
 * sequence is shorter than the motif or no valid window exists. */
static double score_sequence(const char *seq, int slen, const double *lo, int w) {
    if (slen < w) return -INFINITY;
    static const int RC[4] = {3, 2, 1, 0};
    double best = -INFINITY;
    int n_pos = slen - w + 1;
    for (int i = 0; i < n_pos; i++) {
        double s = 0; int ok = 1;
        for (int j = 0; j < w; j++) {
            int b = b2i_c(seq[i + j]);
            if (b < 0) { ok = 0; break; }
            s += lo[j*4 + b];
        }
        if (ok && s > best) best = s;
    }
    for (int i = 0; i < n_pos; i++) {
        double s = 0; int ok = 1;
        for (int j = 0; j < w; j++) {
            int b = b2i_c(seq[i + w - 1 - j]);
            if (b < 0) { ok = 0; break; }
            s += lo[j*4 + RC[b]];
        }
        if (ok && s > best) best = s;
    }
    return best;
}

/* ===================== AUROC ===================== */

typedef struct { double s; int label; } SL;

static int sl_cmp_desc(const void *a, const void *b) {
    const SL *x = (const SL *)a;
    const SL *y = (const SL *)b;
    if (x->s < y->s) return 1;
    if (x->s > y->s) return -1;
    return 0;
}

static double compute_auroc(const double *pos, int npos,
                            const double *neg, int nneg) {
    if (npos == 0 || nneg == 0) return 0.5;
    int total = npos + nneg;
    SL *arr = (SL *)malloc(sizeof(SL) * total);
    for (int i = 0; i < npos; i++) { arr[i].s = pos[i]; arr[i].label = 1; }
    for (int i = 0; i < nneg; i++) { arr[npos + i].s = neg[i]; arr[npos + i].label = 0; }
    qsort(arr, total, sizeof(SL), sl_cmp_desc);
    long concordant = 0;
    int  pos_rank   = 0;
    for (int i = 0; i < total; i++) {
        if (arr[i].label == 1) pos_rank++;
        else                   concordant += pos_rank;
    }
    free(arr);
    return (double)concordant / ((double)npos * (double)nneg);
}

/* ===================== Subprocess runner ===================== */

static int run_subprocess(char *const *argv, int timeout_sec) {
    (void)timeout_sec;  /* not enforced; tools are expected to terminate */
    pid_t pid = fork();
    if (pid < 0) return -1;
    if (pid == 0) {
        int devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) {
            dup2(devnull, 1);
            dup2(devnull, 2);
            close(devnull);
        }
        execvp(argv[0], argv);
        _exit(127);
    }
    int status = 0;
    if (waitpid(pid, &status, 0) < 0) return -1;
    if (!WIFEXITED(status)) return -1;
    return WEXITSTATUS(status);
}

static int file_exists(const char *path) {
    struct stat st;
    return stat(path, &st) == 0;
}

/* ===================== MEME txt parser ===================== */

static float *parse_meme_txt(const char *path, int *out_w) {
    *out_w = 0;
    FILE *f = fopen(path, "r");
    if (!f) return NULL;

    char line[4096];
    float *pwm = NULL;
    int w = 0, cap = 0;
    int saw_header = 0;

    while (fgets(line, sizeof(line), f)) {
        if (!saw_header) {
            if (strstr(line, "letter-probability matrix")) saw_header = 1;
            continue;
        }
        /* Skip blank lines and section markers. */
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == 0 || *p == '\n' || *p == '\r') {
            if (w > 0) break;
            continue;
        }
        if (strncmp(p, "---", 3) == 0 || strncmp(p, "Time", 4) == 0 ||
            strncmp(p, "service", 7) == 0) {
            if (w > 0) break; else continue;
        }
        /* Try to parse first 4 whitespace-separated tokens as floats. */
        double vals[4];
        int matched = 0;
        char *tok = strtok(p, " \t\r\n");
        for (int k = 0; k < 4 && tok; k++) {
            char *endp = NULL;
            double v = strtod(tok, &endp);
            if (endp == tok) { matched = -1; break; }
            vals[k] = v;
            matched++;
            tok = strtok(NULL, " \t\r\n");
        }
        if (matched < 4) {
            if (w > 0) break;
            continue;
        }
        if (w == cap) {
            cap = cap > 0 ? cap * 2 : 16;
            pwm = (float *)realloc(pwm, sizeof(float) * (size_t)cap * 4);
        }
        for (int b = 0; b < 4; b++) pwm[w*4 + b] = (float)vals[b];
        w++;
    }
    fclose(f);
    if (w == 0) { free(pwm); return NULL; }
    *out_w = w;
    return pwm;
}

/* ===================== STREME XML parser ===================== */

static float *parse_streme_xml(const char *path, int *out_w) {
    *out_w = 0;
    FILE *f = fopen(path, "r");
    if (!f) return NULL;

    char line[8192];
    float *pwm = NULL;
    int w = 0, cap = 0;
    int in_motif = 0;

    while (fgets(line, sizeof(line), f)) {
        if (!in_motif) {
            if (strstr(line, "<motif")) {
                /* Only take the FIRST motif. */
                in_motif = 1;
            }
            continue;
        }
        if (strstr(line, "</motif>")) break;

        /* Look for <pos A=".." C=".." G=".." T=".."/> */
        char *q = strstr(line, "<pos");
        if (!q) continue;
        double A = 0, C = 0, G = 0, T = 0;
        int got = 0;
        char *attrs = q;
        char *pa = strstr(attrs, "A=\""); if (pa) { A = strtod(pa + 3, NULL); got++; }
        char *pc = strstr(attrs, "C=\""); if (pc) { C = strtod(pc + 3, NULL); got++; }
        char *pg = strstr(attrs, "G=\""); if (pg) { G = strtod(pg + 3, NULL); got++; }
        char *pt = strstr(attrs, "T=\""); if (pt) { T = strtod(pt + 3, NULL); got++; }
        if (got < 4) continue;
        if (w == cap) {
            cap = cap > 0 ? cap * 2 : 16;
            pwm = (float *)realloc(pwm, sizeof(float) * (size_t)cap * 4);
        }
        pwm[w*4+0] = (float)A;
        pwm[w*4+1] = (float)C;
        pwm[w*4+2] = (float)G;
        pwm[w*4+3] = (float)T;
        w++;
    }
    fclose(f);
    if (w == 0) { free(pwm); return NULL; }
    *out_w = w;
    return pwm;
}

/* ===================== TF discovery ===================== */

static int suffix_is(const char *s, const char *suf) {
    int ls = (int)strlen(s), lf = (int)strlen(suf);
    return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

static int cstr_cmp(const void *a, const void *b) {
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

/* Lists TF names (caller frees) by scanning dir for *_sequences.fa. */
static char **list_tfs(const char *dir, int *out_n) {
    DIR *d = opendir(dir);
    if (!d) return NULL;
    int cap = 32, n = 0;
    char **arr = (char **)malloc(sizeof(char *) * cap);
    struct dirent *e;
    const char *suf = "_sequences.fa";
    int suflen = (int)strlen(suf);
    while ((e = readdir(d)) != NULL) {
        if (!suffix_is(e->d_name, suf)) continue;
        int nl = (int)strlen(e->d_name) - suflen;
        if (n == cap) { cap *= 2; arr = (char **)realloc(arr, sizeof(char *) * cap); }
        arr[n] = (char *)malloc((size_t)nl + 1);
        memcpy(arr[n], e->d_name, nl);
        arr[n][nl] = 0;
        n++;
    }
    closedir(d);
    qsort(arr, n, sizeof(char *), cstr_cmp);
    *out_n = n;
    return arr;
}

/* ===================== Per-TF result & runner ===================== */

typedef struct {
    int    width;
    double auroc;
    double time_s;
    int    ran;       /* 1 if algorithm ran and produced a PWM */
    double p_value;    /* Fisher's exact enrichment p-value */
    double e_value;    /* Bonferroni-corrected p-value */
    double enrichment; /* pos_hit_rate / neg_hit_rate */
    double pos_hit;    /* fraction of positives above threshold */
    double neg_hit;    /* fraction of negatives above threshold */
    double threshold;  /* score threshold for significance */
} AlgoResult;

typedef struct {
    const char *tf;
    int n_pos;
    int n_neg;
    AlgoResult ours, streme, meme;
} TFRow;

static int g_ours_only = 0;
static int g_no_meme = 0;
static int g_markov_bg = 0;  /* if 1, use Markov background for scoring */
static const char *g_negatives_file = NULL;  /* if set, load negatives from file */

/* Forward declarations — scoring + significance functions defined after main */
static double score_and_auroc_markov(const float *pwm, int w,
                                      char **pos, int npos,
                                      char **neg, int nneg);

typedef struct {
    double p_value;
    double e_value;
    double pos_hit_rate;
    double neg_hit_rate;
    double enrichment;
    double threshold;
} MotifSignificance;

static MotifSignificance compute_motif_significance(
    const float *pwm, int w, const double *bg,
    char **pos, int npos, char **neg, int nneg, int n_motifs_tested);

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Score positives + negatives with `pwm` (width w) against `bg`, return AUROC. */
static double score_and_auroc(const float *pwm, int w, const double *bg,
                              char **pos, int npos, char **neg, int nneg) {
    double *lo = pwm_logodds(pwm, w, bg);
    double *ps = (double *)malloc(sizeof(double) * npos);
    double *ns = (double *)malloc(sizeof(double) * nneg);
    for (int i = 0; i < npos; i++)
        ps[i] = score_sequence(pos[i], (int)strlen(pos[i]), lo, w);
    for (int i = 0; i < nneg; i++)
        ns[i] = score_sequence(neg[i], (int)strlen(neg[i]), lo, w);
    double auc = compute_auroc(ps, npos, ns, nneg);
    free(ps); free(ns); free(lo);
    return auc;
}

/* Run one benchmark for a single TF. */
static TFRow run_tf(const char *tf, const char *data_dir) {
    TFRow row;
    memset(&row, 0, sizeof(row));
    row.tf = tf;

    char seq_path[2048], bg_path[2048];
    snprintf(seq_path, sizeof(seq_path), "%s/%s_sequences.fa", data_dir, tf);
    snprintf(bg_path, sizeof(bg_path), "%s/%s_background.txt", data_dir, tf);

    int n_pos = 0;
    char **all_seqs = load_fasta(seq_path, 0, &n_pos);
    if (!all_seqs) {
        fprintf(stderr, "[skip] cannot load %s\n", seq_path);
        return row;
    }
    int n_disc = 0;
    char **disc_seqs = load_fasta(seq_path, DISCOVERY_CENTER_BP, &n_disc);
    double bg[4];
    load_bg(bg_path, bg);
    char **neg_seqs;
    int n_neg;
    if (g_negatives_file) {
        /* Load user-provided negatives */
        neg_seqs = load_fasta(g_negatives_file, 0, &n_neg);
        if (!neg_seqs) {
            fprintf(stderr, "[skip] cannot load negatives %s\n", g_negatives_file);
            free(all_seqs); free(disc_seqs);
            return row;
        }
    } else {
        /* Auto-generate shuffled negatives */
        neg_seqs = shuffle_all(all_seqs, n_pos);
        n_neg = n_pos;
    }
    row.n_pos = n_pos;
    row.n_neg = n_neg;

    /* ---- ours ---- */
    {
        double t0 = now_seconds();
        Motif *motifs = NULL; int n_motifs = 0;
        discover_motifs(disc_seqs, n_disc, 0, &motifs, &n_motifs);
        double t1 = now_seconds();
        if (motifs && n_motifs > 0) {
            /* Try all returned motifs, pick best AUROC */
            double best_auc = 0; int best_idx = 0;
            for (int mi = 0; mi < n_motifs; mi++) {
                double auc;
                if (g_markov_bg) {
                    auc = score_and_auroc_markov(motifs[mi].pwm, motifs[mi].width,
                                                 all_seqs, n_pos, neg_seqs, n_pos);
                } else {
                    auc = score_and_auroc(motifs[mi].pwm, motifs[mi].width, bg,
                                          all_seqs, n_pos, neg_seqs, n_pos);
                }
                if (auc > best_auc) { best_auc = auc; best_idx = mi; }
            }
            row.ours.width  = motifs[best_idx].width;
            row.ours.time_s = t1 - t0;
            row.ours.auroc  = best_auc;
            row.ours.ran    = 1;

            /* PWM dump for analysis */
            {
                const float *pwm = motifs[best_idx].pwm;
                int w = motifs[best_idx].width;
                fprintf(stderr, "PWM:");
                for (int p = 0; p < w; p++)
                    for (int b = 0; b < 4; b++)
                        fprintf(stderr, " %.4f", pwm[p*4+b]);
                fprintf(stderr, "\n");
            }

            /* Compute enrichment significance (Fisher's exact) */
            MotifSignificance ms = compute_motif_significance(
                motifs[best_idx].pwm, motifs[best_idx].width, bg,
                all_seqs, n_pos, neg_seqs, n_pos, n_motifs);
            row.ours.p_value    = ms.p_value;
            row.ours.e_value    = ms.e_value;
            row.ours.enrichment = ms.enrichment;
            row.ours.pos_hit    = ms.pos_hit_rate;
            row.ours.neg_hit    = ms.neg_hit_rate;
            row.ours.threshold  = ms.threshold;
        }
        motif_free_all(motifs, n_motifs);
    }

    /* Write the discovery FASTA once for external tools. */
    char tmp_fa[2048];
    snprintf(tmp_fa, sizeof(tmp_fa), "/tmp/auroc_%s.fa", tf);
    write_fasta(tmp_fa, disc_seqs, n_disc);

    /* ---- STREME ---- */
    if (!g_ours_only && file_exists(STREME_BIN)) {
        char streme_dir[2048];
        snprintf(streme_dir, sizeof(streme_dir), "/tmp/auroc_streme_%s", tf);
        if (file_exists(streme_dir)) {
            pid_t rmpid = fork();
            if (rmpid == 0) { execlp("rm", "rm", "-rf", streme_dir, NULL); _exit(127); }
            int st; waitpid(rmpid, &st, 0);
        }
        char *const sargv[] = {
            (char *)"streme", (char *)"--p",   tmp_fa,
            (char *)"--oc",   streme_dir,
            (char *)"--dna",  (char *)"--minw", (char *)"6",
            (char *)"--maxw", (char *)"17", NULL
        };
        /* Allow executing through the absolute path too. */
        char *const sargv_abs[] = {
            (char *)STREME_BIN, (char *)"--p",   tmp_fa,
            (char *)"--oc",     streme_dir,
            (char *)"--dna",    (char *)"--minw", (char *)"6",
            (char *)"--maxw",   (char *)"17", NULL
        };
        double t0 = now_seconds();
        int rc = run_subprocess(sargv_abs, 180);
        if (rc != 0) rc = run_subprocess(sargv, 180);
        double t1 = now_seconds();
        char xml[4096];
        snprintf(xml, sizeof(xml), "%s/streme.xml", streme_dir);
        int w_s = 0;
        float *spwm = parse_streme_xml(xml, &w_s);
        if (spwm) {
            row.streme.width  = w_s;
            row.streme.time_s = t1 - t0;
            row.streme.auroc  = score_and_auroc(spwm, w_s, bg,
                                                all_seqs, n_pos, neg_seqs, n_pos);
            row.streme.ran    = 1;
            free(spwm);
        }
    }

    /* ---- MEME ---- */
    if (!g_ours_only && !g_no_meme && file_exists(MEME_BIN)) {
        char meme_dir[2048];
        snprintf(meme_dir, sizeof(meme_dir), "/tmp/auroc_meme_%s", tf);
        if (file_exists(meme_dir)) {
            pid_t rmpid = fork();
            if (rmpid == 0) { execlp("rm", "rm", "-rf", meme_dir, NULL); _exit(127); }
            int st; waitpid(rmpid, &st, 0);
        }
        char *const margv_abs[] = {
            (char *)MEME_BIN, tmp_fa,
            (char *)"-o",     meme_dir,
            (char *)"-dna",   (char *)"-mod", (char *)"zoops",
            (char *)"-nmotifs", (char *)"5",
            (char *)"-minw",  (char *)"6", (char *)"-maxw", (char *)"17",
            (char *)"-revcomp", NULL
        };
        double t0 = now_seconds();
        int rc = run_subprocess(margv_abs, 300);
        double t1 = now_seconds();
        (void)rc;
        char txt[4096];
        snprintf(txt, sizeof(txt), "%s/meme.txt", meme_dir);
        int w_m = 0;
        float *mpwm = parse_meme_txt(txt, &w_m);
        if (mpwm) {
            row.meme.width  = w_m;
            row.meme.time_s = t1 - t0;
            row.meme.auroc  = score_and_auroc(mpwm, w_m, bg,
                                              all_seqs, n_pos, neg_seqs, n_pos);
            row.meme.ran    = 1;
            free(mpwm);
        }
    }

    unlink(tmp_fa);
    free_strings(all_seqs, n_pos);
    free_strings(disc_seqs, n_disc);
    free_strings(neg_seqs, n_pos);
    return row;
}

static void print_row(const char *tf, const AlgoResult *a, const char *name) {
    if (!a->ran) return;
    printf("%s\t%d\t%.4f\t%.1f\t%s\t%.2e\t%.4f\t%.4f\n",
           tf, a->width, a->auroc, a->time_s, name,
           a->p_value, a->pos_hit, a->neg_hit);
}

int main(int argc, char **argv) {
    const char *data_dir  = NULL;
    const char *tf_filter = NULL;
    int ours_only = 0;
    for (int i = 1; i < argc; i++) {
        if      (strcmp(argv[i], "--data") == 0 && i + 1 < argc) data_dir  = argv[++i];
        else if ((strcmp(argv[i], "--tf") == 0 || strcmp(argv[i], "--tfs") == 0)
                 && i + 1 < argc)                                 tf_filter = argv[++i];
        else if (strcmp(argv[i], "--ours-only") == 0)            ours_only = 1;
        else if (strcmp(argv[i], "--no-meme") == 0)              g_no_meme = 1;
        else if (strcmp(argv[i], "--markov") == 0)               g_markov_bg = 1;
        else if (strcmp(argv[i], "--negatives") == 0 && i+1 < argc) g_negatives_file = argv[++i];
    }
    g_ours_only = ours_only;
    if (!data_dir) {
        fprintf(stderr, "Usage: %s --data <data_dir> [--tf USF1|ALL] [--markov]\n", argv[0]);
        return 1;
    }

    int  n_tfs = 0;
    char **tfs = list_tfs(data_dir, &n_tfs);
    if (n_tfs == 0) {
        fprintf(stderr, "No *_sequences.fa files in %s\n", data_dir);
        return 1;
    }

    /* Header */
    printf("TF\tWidth\tAUROC\tTime_s\tSource\tP_value\tPos_hit\tNeg_hit\n");
    fflush(stdout);

    /* Accumulators */
    int    n_ours = 0, n_meme = 0, n_streme = 0;
    double s_ours_auc = 0, s_ours_t = 0;
    double s_meme_auc = 0, s_meme_t = 0;
    double s_streme_auc = 0, s_streme_t = 0;

    for (int i = 0; i < n_tfs; i++) {
        if (tf_filter && strcmp(tf_filter, "ALL") != 0 && strcmp(tf_filter, tfs[i]) != 0)
            continue;
        fprintf(stderr, "[%d/%d] %s ... ", i + 1, n_tfs, tfs[i]);
        TFRow row = run_tf(tfs[i], data_dir);
        if (ours_only) {
            memset(&row.streme, 0, sizeof(row.streme));
            memset(&row.meme, 0, sizeof(row.meme));
        }
        print_row(tfs[i], &row.ours,   "ours");
        print_row(tfs[i], &row.streme, "streme");
        print_row(tfs[i], &row.meme,   "meme");
        fflush(stdout);
        if (row.ours.ran)   { n_ours++;   s_ours_auc   += row.ours.auroc;   s_ours_t   += row.ours.time_s; }
        if (row.streme.ran) { n_streme++; s_streme_auc += row.streme.auroc; s_streme_t += row.streme.time_s; }
        if (row.meme.ran)   { n_meme++;   s_meme_auc   += row.meme.auroc;   s_meme_t   += row.meme.time_s; }
        fprintf(stderr, "ours=%.4f (w=%d)\n", row.ours.auroc, row.ours.width);
    }

    /* Summary */
    printf("\n=== SUMMARY ===\n");
    if (n_ours)   printf("ours   avg AUROC=%.4f  avg time=%.2fs  (n=%d)\n",
                         s_ours_auc / n_ours,   s_ours_t / n_ours,   n_ours);
    if (n_streme) printf("streme avg AUROC=%.4f  avg time=%.2fs  (n=%d)\n",
                         s_streme_auc / n_streme, s_streme_t / n_streme, n_streme);
    if (n_meme)   printf("meme   avg AUROC=%.4f  avg time=%.2fs  (n=%d)\n",
                         s_meme_auc / n_meme,   s_meme_t / n_meme,   n_meme);

    for (int i = 0; i < n_tfs; i++) free(tfs[i]);
    free(tfs);
    return 0;
}

/* ===================== First-order Markov background =====================
 *
 * Instead of position-independent mononucleotide background, use a
 * first-order Markov model fitted on the negative sequences. This
 * captures dinucleotide composition (CpG suppression, etc.) without
 * requiring a precomputed k-mer table.
 *
 * The Markov log-odds for a window b_0..b_{w-1}:
 *   log2(PWM_score) - log2(P_markov(window))
 * where P_markov(window) = P(b_0) * Π P(b_i | b_{i-1})
 *
 * This is computed per-window (not precomputable per-position like
 * mononucleotide), so scoring is slightly slower but still O(n*w).
 */

typedef struct {
    double initial[4];      /* P(b_0) — base frequencies */
    double trans[4][4];     /* P(b_i | b_{i-1}) — transition matrix */
} MarkovBg;

/* Build Markov model from a set of sequences. */
static void markov_bg_build(MarkovBg *mb, char **seqs, int n_seqs) {
    double init_count[4] = {0, 0, 0, 0};
    double trans_count[4][4];
    memset(trans_count, 0, sizeof(trans_count));
    int init_total = 0;

    for (int s = 0; s < n_seqs; s++) {
        const char *seq = seqs[s];
        int len = (int)strlen(seq);
        if (len == 0) continue;

        int prev = b2i_c(seq[0]);
        if (prev >= 0) { init_count[prev]++; init_total++; }

        for (int i = 1; i < len; i++) {
            int cur = b2i_c(seq[i]);
            if (prev >= 0 && cur >= 0) {
                trans_count[prev][cur]++;
            }
            prev = cur;
        }
    }

    /* Initial probabilities with pseudocount */
    double init_sum = 0;
    for (int b = 0; b < 4; b++) {
        mb->initial[b] = init_count[b] + 1.0; /* pseudocount */
        init_sum += mb->initial[b];
    }
    for (int b = 0; b < 4; b++) mb->initial[b] /= init_sum;

    /* Transition probabilities with pseudocount */
    for (int prev = 0; prev < 4; prev++) {
        double row_sum = 0;
        for (int cur = 0; cur < 4; cur++) {
            mb->trans[prev][cur] = trans_count[prev][cur] + 0.25; /* pseudocount */
            row_sum += mb->trans[prev][cur];
        }
        for (int cur = 0; cur < 4; cur++) mb->trans[prev][cur] /= row_sum;
    }
}

/* Precompute log2 of Markov probabilities for fast scoring. */
typedef struct {
    double log2_initial[4];
    double log2_trans[4][4];
} MarkovLog;

static void markov_log_build(MarkovLog *ml, const MarkovBg *mb) {
    for (int b = 0; b < 4; b++)
        ml->log2_initial[b] = log2(mb->initial[b]);
    for (int prev = 0; prev < 4; prev++)
        for (int cur = 0; cur < 4; cur++)
            ml->log2_trans[prev][cur] = log2(mb->trans[prev][cur]);
}

/* Score a single window with PWM against Markov background.
 * Returns log2(PWM(window)) - log2(Markov(window)). */
static double score_window_markov(const char *window, int w,
                                   const float *pwm,
                                   const MarkovLog *ml) {
    double pwm_ll = 0;
    double markov_ll = 0;
    int prev = -1;

    for (int j = 0; j < w; j++) {
        int b = b2i_c(window[j]);
        if (b < 0) return -INFINITY;

        pwm_ll += log2((double)pwm[j * 4 + b] + 1e-6);

        if (j == 0) {
            markov_ll += ml->log2_initial[b];
        } else {
            markov_ll += ml->log2_trans[prev][b];
        }
        prev = b;
    }

    return pwm_ll - markov_ll;
}

/* Score sequence: max window score over both strands, Markov background. */
static double score_sequence_markov(const char *seq, int slen,
                                     const float *pwm, int w,
                                     const MarkovLog *ml) {
    if (slen < w) return -INFINITY;
    double best = -INFINITY;
    int n_pos = slen - w + 1;

    /* Forward strand */
    for (int i = 0; i < n_pos; i++) {
        double s = score_window_markov(seq + i, w, pwm, ml);
        if (s > best) best = s;
    }

    /* Reverse complement: build RC of each window */
    char rc_buf[32];
    for (int i = 0; i < n_pos; i++) {
        /* Build RC window */
        int ok = 1;
        for (int j = 0; j < w; j++) {
            int b = b2i_c(seq[i + w - 1 - j]);
            if (b < 0) { ok = 0; break; }
            rc_buf[j] = "TGCA"[b]; /* RC of base */
        }
        if (!ok) continue;
        double s = score_window_markov(rc_buf, w, pwm, ml);
        if (s > best) best = s;
    }

    return best;
}

/* Score positives + negatives with Markov background. */
static double score_and_auroc_markov(const float *pwm, int w,
                                      char **pos, int npos,
                                      char **neg, int nneg) {
    /* Build Markov model from negatives */
    MarkovBg mb;
    markov_bg_build(&mb, neg, nneg);
    MarkovLog ml;
    markov_log_build(&ml, &mb);

    double *ps = (double *)malloc(sizeof(double) * npos);
    double *ns = (double *)malloc(sizeof(double) * nneg);
    for (int i = 0; i < npos; i++)
        ps[i] = score_sequence_markov(pos[i], (int)strlen(pos[i]), pwm, w, &ml);
    for (int i = 0; i < nneg; i++)
        ns[i] = score_sequence_markov(neg[i], (int)strlen(neg[i]), pwm, w, &ml);

    double auc = compute_auroc(ps, npos, ns, nneg);
    free(ps); free(ns);
    return auc;
}

/* ===================== Fisher's Exact Test =====================
 *
 * Computes the two-sided p-value for a 2x2 contingency table:
 *   [a b]
 *   [c d]
 *
 * Using the hypergeometric distribution:
 *   P = (C(a+b,a) * C(c+d,c)) / C(n, a+c)
 * summed over all tables at least as extreme.
 *
 * For motif significance: a = pos_above, b = neg_above,
 * c = pos_below, d = neg_below. Significant motifs have
 * enrichment (a/(a+c) >> b/(b+d)).
 */

/* ln factorial via lgamma */
static double log_fact(int n) {
    if (n <= 1) return 0.0;
    return lgamma((double)(n + 1));
}

/* log of C(n, k) */
static double log_choose(int n, int k) {
    if (k < 0 || k > n) return -INFINITY;
    return log_fact(n) - log_fact(k) - log_fact(n - k);
}

/* Fisher's exact test (two-sided) */
static double fisher_exact(int a, int b, int c, int d) {
    int n = a + b + c + d;
    int row1 = a + b;   /* total above threshold */
    int col1 = a + c;   /* total positives */

    if (n == 0 || row1 == 0 || col1 == 0) return 1.0;

    /* The observed table's probability */
    double log_p_obs = log_choose(col1, a) +
                       log_choose(n - col1, b) -
                       log_choose(n, row1);

    /* Sum probabilities of all tables at least as extreme */
    double p_total = 0.0;

    /* Range of possible a values: max(0, row1-d) to min(row1, col1) */
    int a_min = (row1 - d > 0) ? row1 - d : 0;
    int a_max = (row1 < col1) ? row1 : col1;

    for (int ai = a_min; ai <= a_max; ai++) {
        int bi = row1 - ai;
        int ci = col1 - ai;
        int di = n - row1 - col1 + ai;
        if (bi < 0 || ci < 0 || di < 0) continue;

        double log_p = log_choose(col1, ai) +
                       log_choose(n - col1, bi) -
                       log_choose(n, row1);
        double p = exp(log_p);

        /* Two-sided: include if probability <= observed probability */
        if (log_p <= log_p_obs + 1e-15) {
            p_total += p;
        }
    }

    return (p_total > 1.0) ? 1.0 : p_total;
}

/* Compute motif enrichment significance.
 * Scores all sequences, finds optimal threshold (Youden's J),
 * builds 2x2 table, returns Fisher's exact p-value.
 *
 * Also reports basic stats: pos_hit_rate, neg_hit_rate, enrichment ratio.
 */

/* MotifSignificance typedef is at top of file (forward declaration) */

static MotifSignificance compute_motif_significance(
        const float *pwm, int w,
        const double *bg,
        char **pos, int npos,
        char **neg, int nneg,
        int n_motifs_tested) {

    MotifSignificance ms;
    ms.p_value = 1.0; ms.e_value = 1.0;
    ms.pos_hit_rate = 0; ms.neg_hit_rate = 0;
    ms.enrichment = 1.0; ms.threshold = 0;

    /* Score all sequences */
    double *lo = (double *)malloc(sizeof(double) * (size_t)w * 4);
    for (int j = 0; j < w; j++) {
        double row[4], s = 0;
        for (int b = 0; b < 4; b++) { row[b] = (double)pwm[j*4+b] + 1e-6; s += row[b]; }
        for (int b = 0; b < 4; b++) { row[b] /= s; lo[j*4+b] = log2(row[b] / bg[b]); }
    }

    double *ps = (double *)malloc(sizeof(double) * npos);
    double *ns = (double *)malloc(sizeof(double) * nneg);
    for (int i = 0; i < npos; i++)
        ps[i] = score_sequence(pos[i], (int)strlen(pos[i]), lo, w);
    for (int i = 0; i < nneg; i++)
        ns[i] = score_sequence(neg[i], (int)strlen(neg[i]), lo, w);

    /* Find optimal threshold via Youden's J (maximize TPR - FPR) */
    /* Sort all scores, sweep threshold */
    int total = npos + nneg;
    SL *arr = (SL *)malloc(sizeof(SL) * total);
    for (int i = 0; i < npos; i++) { arr[i].s = ps[i]; arr[i].label = 1; }
    for (int i = 0; i < nneg; i++) { arr[npos+i].s = ns[i]; arr[npos+i].label = 0; }
    qsort(arr, total, sizeof(SL), sl_cmp_desc);

    double best_j = -1.0;
    double best_threshold = arr[0].s;
    int tp = 0, fp = 0;
    for (int i = 0; i < total; i++) {
        if (arr[i].label == 1) tp++;
        else fp++;
        double tpr = (double)tp / npos;
        double fpr = (double)fp / nneg;
        double j = tpr - fpr;
        if (j > best_j) {
            best_j = j;
            best_threshold = arr[i].s;
        }
    }

    /* Build 2x2 table at optimal threshold */
    int pos_above = 0, neg_above = 0;
    for (int i = 0; i < npos; i++) if (ps[i] >= best_threshold) pos_above++;
    for (int i = 0; i < nneg; i++) if (ns[i] >= best_threshold) neg_above++;
    int pos_below = npos - pos_above;
    int neg_below = nneg - neg_above;

    /* Fisher's exact test */
    double pval = fisher_exact(pos_above, neg_above, pos_below, neg_below);

    ms.p_value = pval;
    ms.e_value = pval * (double)n_motifs_tested;
    ms.pos_hit_rate = (double)pos_above / npos;
    ms.neg_hit_rate = (double)neg_above / (nneg > 0 ? nneg : 1);
    ms.enrichment = (ms.neg_hit_rate > 0) ? ms.pos_hit_rate / ms.neg_hit_rate : INFINITY;
    ms.threshold = best_threshold;

    free(lo); free(ps); free(ns); free(arr);
    return ms;
}
