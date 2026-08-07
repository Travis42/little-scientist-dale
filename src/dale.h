#ifndef MOTIF_H
#define MOTIF_H

#ifdef __cplusplus
extern "C" {
#endif

#define MOTIF_MAX_RETURN 5
#define MOTIF_INFO_LEN   128

/* A discovered motif. PWM is row-major: pwm[i*4 + b] for position i, base b
 * where b = 0 (A), 1 (C), 2 (G), 3 (T). */
typedef struct {
    float *pwm;             /* width * 4 floats, row-major (A,C,G,T per row) */
    int    width;           /* number of positions */
    float  ic;              /* mean information content (bits) against bg */
    char   info[MOTIF_INFO_LEN]; /* human-readable info string */
} Motif;

/* Discover motifs in a set of DNA sequences (upper-case A/C/G/T strings).
 *
 *   sequences   : array of NUL-terminated DNA strings
 *   n_sequences : number of sequences
 *   k_hint      : if >0, restrict to that single width; if 0 use full 6..17 range
 *   out_motifs  : receives malloc'd array of Motif (caller frees each
 *                 motif's `pwm` and the array itself with motif_free_all)
 *   out_n_motifs: receives number of motifs (<= MOTIF_MAX_RETURN)
 *
 * Returns 0 on success, non-zero on error. */
int discover_motifs(char **sequences, int n_sequences, int k_hint,
                    Motif **out_motifs, int *out_n_motifs);

/* Convenience: free an array of Motif returned by discover_motifs. */
void motif_free_all(Motif *motifs, int n_motifs);

#ifdef __cplusplus
}
#endif

#endif /* MOTIF_H */
