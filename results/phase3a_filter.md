## Phase 3A, item 1: cheaper seed selection and filter_auto

### A. Seed-selection variants (E. coli 536, m = 100, pure filter, mean over reads)

Selection time is the mean time spent choosing seeds per query (100 reads, best of 3 runs); hits are the mean seed hits of the chosen seeds, so a variant that is cheaper but chooses worse seeds shows up as more hits.

| k | original DP (all lengths 1..11, per-seed lookups) | prefix tables, all lengths 1..11 | tables + length cap 3 | tables + length cap 2 (default) | tables + length cap 1 |
|---|---|---|---|---|---|
| 1 | 17.3 us, 1 hits, 0.020 ms | 8.6 us, 1 hits, 0.010 ms | 4.5 us, 1 hits, 0.006 ms | 4.1 us, 1 hits, 0.006 ms | 2.9 us, 1 hits, 0.005 ms |
| 4 | 20.2 us, 3 hits, 0.024 ms | 11.2 us, 3 hits, 0.015 ms | 5.8 us, 3 hits, 0.009 ms | 4.7 us, 3 hits, 0.008 ms | 3.6 us, 3 hits, 0.007 ms |
| 8 | 26.6 us, 16 hits, 0.049 ms | 15.2 us, 16 hits, 0.039 ms | 6.8 us, 16 hits, 0.029 ms | 6.8 us, 16 hits, 0.032 ms | 4.7 us, 16 hits, 0.029 ms |
| 10 | 34.4 us, 215 hits, 0.403 ms | 19.3 us, 215 hits, 0.375 ms | 11.7 us, 215 hits, 0.355 ms | 10.0 us, 215 hits, 0.366 ms | 7.4 us, 222 hits, 0.384 ms |
| 12 | 49.3 us, 1760 hits, 2.997 ms | 39.3 us, 1760 hits, 3.511 ms | 19.2 us, 1760 hits, 3.073 ms | 16.2 us, 1760 hits, 3.149 ms | 9.9 us, 1801 hits, 3.100 ms |
| 15 | 73.6 us, 14685 hits, 22.648 ms | 46.6 us, 14685 hits, 22.928 ms | 27.5 us, 14685 hits, 21.580 ms | 25.6 us, 14685 hits, 22.874 ms | 18.0 us, 14770 hits, 22.658 ms |

### B. filter_auto versus filter and filter_fa (fallback active; mean ms per query, best of 3 runs)

`auto/best` is the mean time of filter_auto divided by the smaller of the means of filter and filter_fa (target: at most 1.05). Choices are counts of queries per decision: skip = DP not run (a subset of eq), eq, fa, scan.

| text | setting | filter | filter_fa | filter_auto | auto/best | choices (skip / eq / fa / scan) |
|---|---|---|---|---|---|---|
| E. coli 536 | m=100, k=1 | 0.0075 | 0.0061 | 0.0053 | 0.87 | 8 / 9 / 91 / 0 |
| E. coli 536 | m=100, k=2 | 0.0136 | 0.0059 | 0.0058 | 0.99 | 2 / 2 / 98 / 0 |
| E. coli 536 | m=100, k=4 | 0.0220 | 0.0078 | 0.0081 | 1.05 | 0 / 0 / 100 / 0 |
| E. coli 536 | m=100, k=6 | 0.0278 | 0.0110 | 0.0112 | 1.02 | 0 / 0 / 100 / 0 |
| E. coli 536 | m=100, k=8 | 0.0508 | 0.0302 | 0.0291 | 0.96 | 0 / 7 / 93 / 0 |
| E. coli 536 | m=100, k=10 | 0.8039 | 0.3715 | 0.3997 | 1.08 | 0 / 1 / 99 / 0 |
| E. coli 536 | m=100, k=11 | 2.2446 | 1.2523 | 1.2484 | 1.00 | 0 / 1 / 99 / 0 |
| E. coli 536 | m=100, k=12 | 7.8155 | 3.5026 | 3.0968 | 0.88 | 0 / 0 / 100 / 0 |
| E. coli 536 | m=100, k=13 | 12.0550 | 6.7190 | 6.5412 | 0.97 | 0 / 1 / 99 / 0 |
| E. coli 536 | m=100, k=14 | 28.2478 | 13.3296 | 12.8947 | 0.97 | 0 / 0 / 100 / 0 |
| E. coli 536 | m=100, k=15 | 35.8575 | 22.0568 | 23.0967 | 1.05 | 0 / 0 / 100 / 0 |
| E. coli 536 | m=100, k=17 | 66.3081 | 48.0851 | 49.4433 | 1.03 | 0 / 0 / 100 / 0 |
| E. coli 536 | m=100, k=20 | 66.4618 | 67.1910 | 64.9385 | 0.98 | 0 / 0 / 2 / 98 |
| E. coli 536 | m=200, k=10 | 0.1520 | 0.0256 | 0.0242 | 0.95 | 0 / 0 / 100 / 0 |
| E. coli 536 | m=400, k=20 | 1.0157 | 0.0923 | 0.0919 | 1.00 | 0 / 0 / 100 / 0 |
| E. coli 536 | m=800, k=40 | 7.2595 | 0.3951 | 0.4208 | 1.07 | 0 / 0 / 100 / 0 |
| synthetic | m=100, k=2 | 0.0082 | 0.0067 | 0.0056 | 0.83 | 1 / 1 / 99 / 0 |
| synthetic | m=100, k=8 | 0.0306 | 0.0204 | 0.0207 | 1.02 | 0 / 7 / 93 / 0 |
| synthetic | m=100, k=10 | 0.5589 | 0.3446 | 0.3289 | 0.95 | 0 / 16 / 84 / 0 |
| synthetic | m=100, k=12 | 5.8112 | 2.8739 | 2.8752 | 1.00 | 0 / 1 / 99 / 0 |
| synthetic | m=100, k=15 | 29.3480 | 22.6892 | 22.0704 | 0.97 | 0 / 3 / 97 / 0 |
| stress | low-complexity edge, k=6 | 0.4552 | 0.0381 | 0.0392 | 1.03 | 0 / 1 / 99 / 0 |
| stress | low-complexity inside, k=6 | 0.7407 | 0.5956 | 0.5677 | 0.95 | 0 / 0 / 100 / 0 |
| stress | clustered, k=6 | 0.0425 | 0.0208 | 0.0197 | 0.94 | 0 / 0 / 100 / 0 |
| stress | low-complexity edge, k=10 | 1.3373 | 2.2465 | 1.2881 | 0.96 | 0 / 49 / 51 / 0 |
| stress | low-complexity inside, k=10 | 1.7407 | 1.4518 | 1.3958 | 0.96 | 0 / 17 / 83 / 0 |
| stress | clustered, k=10 | 1.0221 | 0.4524 | 0.4539 | 1.00 | 0 / 0 / 100 / 0 |
| stress | low-complexity edge, k=15 | 31.9837 | 22.2839 | 21.9175 | 0.98 | 0 / 14 / 86 / 0 |
| stress | low-complexity inside, k=15 | 24.9067 | 18.6762 | 18.2589 | 0.98 | 0 / 14 / 86 / 0 |
| stress | clustered, k=15 | 37.6503 | 23.3706 | 23.1225 | 0.99 | 0 / 1 / 99 / 0 |

Worst auto/best ratio: 1.08 (E. coli 536, m=100, k=10); 2 of 30 settings exceed 1.05: E. coli 536 m=100, k=10 (1.08); E. coli 536 m=800, k=40 (1.07).

### Run-to-run noise check (E. coli 536, 300 reads, 6 interleaved runs)

| setting | filter_fa runs: min to max (ms) | filter_auto runs: min to max (ms) | spread of the same algorithm | ratio of medians auto / fa | ratio of minima |
|---|---|---|---|---|---|
| m=100, k=10 | 0.408 to 0.580 | 0.410 to 0.596 | 42% | 0.887 | 1.005 |
| m=800, k=40 | 0.407 to 0.659 | 0.443 to 0.651 | 62% | 1.296 | 1.089 |

The same algorithm varies by more than 5% between identical runs, so single-run ratios of 1.05 to 1.08 are within noise.
