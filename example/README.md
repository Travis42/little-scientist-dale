# Example Dataset

11 transcription factor ChIP-seq profiles from ENCODE K562.

Each file contains 200 peak-centered 500bp sequences in FASTA format.

**Source:** ENCODE Project Consortium (open access). These are standard ENCODE
narrowPeak datasets from K562 cell line ChIP-seq experiments.

**TFs included:**
- ATF3, CTCF, ELF1, GABPA, GATA2, JUN, MAX, MYC, RXRA, SP1, TAL1

These 11 TFs show a range of results: DALE significantly outperforms
STREME on several (RXRA, GABPA, CTCF, SP1), closely matches on others (GATA2, JUN),
and STREME wins on one (ELF1). This was chosen to be an honest representation.

## Quick test

```
./motif --data example/ --markov
```

Expected runtime: ~15 seconds.
