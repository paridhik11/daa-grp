## Phase 3A, item 1: cheaper seed selection and filter_auto

### A. Seed-selection variants (E. coli 536, m = 100, pure filter, mean over reads)

Selection time is the mean time spent choosing seeds per query (100 reads, best of 3 runs); hits are the mean seed hits of the chosen seeds, so a variant that is cheaper but chooses worse seeds shows up as more hits.

| k | original DP (all lengths 1..11, per-seed lookups) | prefix tables, all lengths 1..11 | tables + length cap 3 | tables + length cap 2 (default) | tables + length cap 1 |
|---|---|---|---|---|---|
| 1 | 15.2 us, 1 hits, 0.018 ms | 9.3 us, 1 hits, 0.012 ms | 5.5 us, 1 hits, 0.007 ms | 4.2 us, 1 hits, 0.006 ms | 3.3 us, 1 hits, 0.005 ms |
| 4 | 22.2 us, 3 hits, 0.026 ms | 12.6 us, 3 hits, 0.017 ms | 6.7 us, 3 hits, 0.010 ms | 5.2 us, 3 hits, 0.009 ms | 4.2 us, 3 hits, 0.008 ms |
| 8 | 30.3 us, 16 hits, 0.057 ms | 17.6 us, 16 hits, 0.045 ms | 9.4 us, 16 hits, 0.038 ms | 6.9 us, 16 hits, 0.033 ms | 5.4 us, 16 hits, 0.033 ms |
| 10 | 38.1 us, 215 hits, 0.447 ms | 22.7 us, 215 hits, 0.434 ms | 14.0 us, 215 hits, 0.424 ms | 11.6 us, 215 hits, 0.423 ms | 9.0 us, 222 hits, 0.430 ms |
| 12 | 56.3 us, 1760 hits, 3.435 ms | 31.4 us, 1760 hits, 3.422 ms | 20.7 us, 1760 hits, 3.439 ms | 15.9 us, 1760 hits, 3.357 ms | 10.9 us, 1801 hits, 3.531 ms |
| 15 | 83.1 us, 14685 hits, 24.628 ms | 50.1 us, 14685 hits, 24.708 ms | 33.2 us, 14685 hits, 24.620 ms | 28.0 us, 14685 hits, 24.790 ms | 20.5 us, 14770 hits, 25.076 ms |

### B. filter_auto versus filter and filter_fa (fallback active; mean ms per query, best of 3 runs)

`auto/best` is the mean time of filter_auto divided by the smaller of the means of filter and filter_fa (target: at most 1.05). Choices are counts of queries per decision: skip = DP not run (a subset of eq), eq, fa, scan.

| text | setting | filter | filter_fa | filter_auto | auto/best | choices (skip / eq / fa / scan) |
|---|---|---|---|---|---|---|
| E. coli 536 | m=100, k=1 | 0.0091 | 0.0062 | 0.0057 | 0.91 | 12 / 12 / 87 / 0 |
| E. coli 536 | m=100, k=2 | 0.0145 | 0.0070 | 0.0070 | 0.99 | 2 / 2 / 98 / 0 |
| E. coli 536 | m=100, k=4 | 0.0257 | 0.0092 | 0.0091 | 0.99 | 0 / 0 / 98 / 0 |
| E. coli 536 | m=100, k=6 | 0.0326 | 0.0125 | 0.0125 | 1.00 | 0 / 0 / 72 / 0 |
| E. coli 536 | m=100, k=8 | 0.0586 | 0.0358 | 0.0257 | 0.72 | 0 / 1 / 1 / 0 |
| E. coli 536 | m=100, k=10 | 0.9166 | 0.4260 | 0.0663 | 0.16 | 0 / 0 / 0 / 0 |
| E. coli 536 | m=100, k=11 | 2.5446 | 1.3770 | 0.1417 | 0.10 | 0 / 0 / 0 / 0 |
| E. coli 536 | m=100, k=12 | 8.3791 | 3.3939 | 0.3446 | 0.10 | 0 / 0 / 0 / 0 |
| E. coli 536 | m=100, k=13 | 13.3190 | 7.4962 | 0.8870 | 0.12 | 0 / 0 / 0 / 0 |
| E. coli 536 | m=100, k=14 | 31.0514 | 13.9657 | 2.0552 | 0.15 | 0 / 0 / 0 / 0 |
| E. coli 536 | m=100, k=15 | 38.7411 | 24.3493 | 3.8209 | 0.16 | 0 / 0 / 0 / 0 |
| E. coli 536 | m=100, k=17 | 70.8597 | 53.1276 | 16.6739 | 0.31 | 0 / 0 / 0 / 0 |
| E. coli 536 | m=100, k=20 | 73.2853 | 72.2102 | 72.4674 | 1.00 | 0 / 0 / 0 / 93 |
| E. coli 536 | m=200, k=10 | 0.1637 | 0.0276 | 0.0374 | 1.36 | 0 / 0 / 34 / 0 |
| E. coli 536 | m=400, k=20 | 1.1163 | 0.1043 | 0.1371 | 1.32 | 0 / 0 / 4 / 0 |
| E. coli 536 | m=800, k=40 | 7.7142 | 0.4405 | 0.4989 | 1.13 | 0 / 0 / 0 / 0 |
| synthetic | m=100, k=2 | 0.0078 | 0.0070 | 0.0065 | 0.93 | 2 / 2 / 98 / 0 |
| synthetic | m=100, k=8 | 0.0345 | 0.0241 | 0.0238 | 0.99 | 0 / 0 / 6 / 0 |
| synthetic | m=100, k=10 | 0.6410 | 0.3748 | 0.0601 | 0.16 | 0 / 0 / 0 / 0 |
| synthetic | m=100, k=12 | 6.2198 | 3.1904 | 0.3022 | 0.09 | 0 / 0 / 0 / 0 |
| synthetic | m=100, k=15 | 32.8081 | 24.7178 | 4.3550 | 0.18 | 0 / 0 / 0 / 0 |
| stress | low-complexity edge, k=6 | 0.4786 | 0.0415 | 0.0388 | 0.94 | 0 / 0 / 41 / 0 |
| stress | low-complexity inside, k=6 | 0.7732 | 0.5650 | 0.5982 | 1.06 | 0 / 0 / 95 / 0 |
| stress | clustered, k=6 | 0.0455 | 0.0226 | 0.0263 | 1.16 | 0 / 0 / 13 / 0 |
| stress | low-complexity edge, k=10 | 1.4636 | 2.4636 | 0.4377 | 0.30 | 0 / 0 / 0 / 0 |
| stress | low-complexity inside, k=10 | 1.7684 | 1.5622 | 1.1296 | 0.72 | 0 / 1 / 52 / 0 |
| stress | clustered, k=10 | 1.0861 | 0.4905 | 0.0919 | 0.19 | 0 / 0 / 1 / 0 |
| stress | low-complexity edge, k=15 | 33.1186 | 24.8621 | 6.5758 | 0.26 | 0 / 0 / 0 / 0 |
| stress | low-complexity inside, k=15 | 26.1813 | 19.3731 | 9.6869 | 0.50 | 0 / 0 / 0 / 0 |
| stress | clustered, k=15 | 37.5153 | 25.2135 | 4.2124 | 0.17 | 0 / 1 / 1 / 0 |

Worst auto/best ratio: 1.36 (E. coli 536, m=200, k=10); 5 of 30 settings exceed 1.05: E. coli 536 m=200, k=10 (1.36); E. coli 536 m=400, k=20 (1.32); stress clustered, k=6 (1.16); E. coli 536 m=800, k=40 (1.13); stress low-complexity inside, k=6 (1.06).

### Run-to-run noise check (E. coli 536, 300 reads, 6 interleaved runs)

| setting | filter_fa runs: min to max (ms) | filter_auto runs: min to max (ms) | spread of the same algorithm | ratio of medians auto / fa | ratio of minima |
|---|---|---|---|---|---|
| m=100, k=10 | 0.441 to 0.507 | 0.067 to 0.070 | 15% | 0.154 | 0.152 |
| m=800, k=40 | 0.450 to 0.529 | 0.506 to 0.539 | 18% | 1.130 | 1.124 |

The same algorithm varies by more than 5% between identical runs, so single-run ratios of 1.05 to 1.08 are within noise.
