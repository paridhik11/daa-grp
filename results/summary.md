# Benchmark summary

Main text: **E. coli 536** (see `data/SOURCES.txt`; no K-12 genome was reachable). Search time is per pattern, single thread; index build is reported separately.

Edlib and WFA2 return the BEST match per read (minimum distance, ends at that minimum), not all matches within k, and include Python call overhead; compare them on time only.

Cost-model constants (ns), calibrated on synthetic genome_synth.txt and evaluated on E. coli 536: c_hit=82.1, c_win=404, c_bp=5.10, a_scan=5.21, c_sel=45.2


Headline experiment (E. coli 536, m = 100, 100 reads per setting): values are the **median [10th, 90th percentile] over reads** of per-read query time; other experiments use means over 20 reads (Sellers 1, WFA2 5).


## Frequency-aware vs equal split (pure filter, no fallback), E. coli 536

| k | hits equal | hits fa | hit ratio (medians) | bp verified equal | bp verified fa | time equal (ms) | time fa (ms) | fa speedup vs equal | fa speedup vs Myers |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 5 [2, 10] | 1 [1, 2] | 5.0x | 410 [104, 927] | 102 [101, 103] | 0.006 [0.002, 0.015] | 0.022 [0.020, 0.025] | 0.29x | 2410x |
| 2 | 9 [5, 16] | 1 [1, 2] | 9.0x | 736 [404, 1477] | 104 [103, 208] | 0.011 [0.006, 0.022] | 0.027 [0.023, 0.033] | 0.42x | 1977x |
| 3 | 12 [6, 18] | 2 [1, 3] | 6.0x | 1041 [526, 1739] | 107 [105, 214] | 0.016 [0.008, 0.027] | 0.028 [0.026, 0.032] | 0.55x | 1840x |
| 4 | 14 [9, 25] | 3 [1, 4] | 4.5x | 1242 [641, 2424] | 110 [107, 322] | 0.018 [0.010, 0.034] | 0.036 [0.027, 0.043] | 0.50x | 1458x |
| 6 | 18 [11, 31] | 4 [3, 7] | 4.5x | 1632 [1008, 3145] | 237 [114, 666] | 0.024 [0.015, 0.044] | 0.032 [0.028, 0.038] | 0.76x | 1656x |
| 8 | 24 [14, 68] | 15 [8, 25] | 1.6x | 2437 [1381, 6990] | 1400 [708, 2439] | 0.033 [0.018, 0.098] | 0.061 [0.046, 0.086] | 0.55x | 863x |
| 10 | 326 [216, 1184] | 212 [115, 321] | 1.5x | 38940 [26171, 134831] | 24740 [13387, 37104] | 0.458 [0.311, 1.646] | 0.328 [0.198, 0.479] | 1.40x | 159x |
| 11 | 1260 [978, 1592] | 667 [347, 1068] | 1.9x | 150050 [116381, 187210] | 80439 [43311, 125188] | 1.712 [1.313, 2.164] | 0.954 [0.539, 1.493] | 1.79x | 55x |
| 12 | 4828 [1813, 6220] | 1656 [1071, 2639] | 2.9x | 558434 [227061, 711286] | 204320 [133726, 310830] | 6.470 [2.573, 8.315] | 2.408 [1.562, 3.705] | 2.69x | 22x |
| 13 | 5494 [4530, 19394] | 3774 [2303, 6031] | 1.5x | 654841 [544943, 1873433] | 456108 [288378, 689800] | 7.362 [6.173, 20.888] | 5.438 [3.586, 8.642] | 1.35x | 10x |
| 14 | 20796 [13438, 26218] | 7370 [5305, 11344] | 2.8x | 2041166 [1413272, 2405927] | 867459 [639550, 1244160] | 23.601 [16.231, 27.931] | 10.525 [7.708, 15.013] | 2.24x | 5x |
| 15 | 23647 [19932, 28911] | 14000 [9854, 20383] | 1.7x | 2277172 [2022391, 2639100] | 1523332 [1110920, 1984815] | 27.164 [23.990, 35.714] | 18.545 [13.984, 27.494] | 1.46x | 3x |
| 17 | 97332 [83921, 111570] | 41428 [28912, 59795] | 2.3x | 4570596 [4396164, 4663902] | 3291704 [2696102, 3903638] | 54.041 [50.739, 59.048] | 39.579 [29.387, 47.390] | 1.37x | 1x |
| 20 | 435776 [128342, 492056] | 144428 [108145, 209961] | 3.0x | 4938868 [4816463, 4938916] | 4849136 [4684707, 4914649] | 82.866 [68.150, 102.836] | 59.340 [55.401, 64.396] | 1.40x | 1x |
| 30 | 2504724 [2330961, 2638079] | 1747860 [1409690, 2349176] | 1.4x | 4938920 [4938920, 4938920] | 4938920 [4938920, 4938920] | 275.081 [257.218, 300.203] | 202.880 [166.620, 270.323] | 1.36x | 0x |

## Cost model (predicted vs measured), E. coli 536

Fit: NNLS (scipy.optimize.nnls) on synthetic genome_synth.txt; evaluated here on E. coli 536 (different text). c_win = 404 ns (unconstrained least squares gave 404, so the sign constraint was not active). Hit counts come from the index; windows and merged bp are estimated assuming uniform hit placement.

- Crossover filter vs Myers: measured k/m = 0.169 (k = 16.9), predicted k/m = 0.166 (k = 16.6)
- Crossover filter_fa vs Myers: measured k/m = 0.191 (k = 19.1), predicted k/m = 0.184 (k = 18.4)
- filter time prediction: median factor error 1.21x, worst 1.24x over 14 values of k (measured median > 0.01 ms; prediction ratio uses medians)
- filter_fa time prediction: median factor error 1.24x, worst 2.38x over 15 values of k (measured median > 0.01 ms; prediction ratio uses medians)
- Myers scan prediction: median factor error 1.02x, worst 1.04x

| k | filter measured (ms) | predicted | filter_fa measured | predicted | Myers measured | predicted |
|---|---|---|---|---|---|---|
| 1 | 0.006 [0.002, 0.015] | 0.008 | 0.022 [0.020, 0.025] | 0.051 | 51.9 [49.9, 53.7] | 51.5 |
| 2 | 0.011 [0.006, 0.022] | 0.014 | 0.027 [0.023, 0.033] | 0.052 | 52.4 [49.9, 57.3] | 51.5 |
| 3 | 0.016 [0.008, 0.027] | 0.019 | 0.028 [0.026, 0.032] | 0.052 | 52.3 [50.4, 55.1] | 51.5 |
| 4 | 0.018 [0.010, 0.034] | 0.021 | 0.036 [0.027, 0.043] | 0.054 | 51.8 [49.8, 54.0] | 51.5 |
| 6 | 0.024 [0.015, 0.044] | 0.030 | 0.032 [0.028, 0.038] | 0.057 | 52.4 [50.8, 58.0] | 51.5 |
| 8 | 0.033 [0.018, 0.098] | 0.040 | 0.061 [0.046, 0.086] | 0.075 | 52.3 [50.1, 56.0] | 51.5 |
| 10 | 0.458 [0.311, 1.646] | 0.559 | 0.328 [0.198, 0.479] | 0.410 | 52.1 [50.2, 54.0] | 51.5 |
| 11 | 1.712 [1.313, 2.164] | 2.127 | 0.954 [0.539, 1.493] | 1.193 | 52.3 [50.9, 53.3] | 51.5 |
| 12 | 6.470 [2.573, 8.315] | 7.835 | 2.408 [1.562, 3.705] | 2.920 | 52.9 [50.7, 58.5] | 51.5 |
| 13 | 7.362 [6.173, 20.888] | 9.073 | 5.438 [3.586, 8.642] | 6.419 | 51.9 [50.3, 53.6] | 51.5 |
| 14 | 23.601 [16.231, 27.931] | 27.619 | 10.525 [7.708, 15.013] | 11.919 | 53.1 [50.7, 61.7] | 51.5 |
| 15 | 27.164 [23.990, 35.714] | 30.562 | 18.545 [13.984, 27.494] | 20.819 | 53.4 [50.7, 59.5] | 51.5 |
| 17 | 54.041 [50.739, 59.048] | 57.621 | 39.579 [29.387, 47.390] | 43.042 | 52.8 [50.6, 56.8] | 51.5 |
| 20 | 82.866 [68.150, 102.836] | 86.172 | 59.340 [55.401, 64.396] | 62.425 | 52.5 [50.4, 57.0] | 51.5 |
| 30 | 275.081 [257.218, 300.203] | 256.099 | 202.880 [166.620, 270.323] | 193.985 | 52.2 [50.1, 54.2] | 51.5 |

## Index build vs query time and break-even, E. coli 536

Q* = index_ms / (Myers median ms per query - filter median ms per query); the number of queries after which building the index plus filtering beats repeated plain Myers scans (never, if the filter is not faster). Query times are medians over 100 reads (with the cost-model fallback active); the index build time is one measurement per setting.

| k | index build (ms) | Myers (ms/query) | filter (ms/query) | Q* filter | filter_fa (ms/query) | Q* filter_fa |
|---|---|---|---|---|---|---|
| 1 | 112 | 51.9 | 0.007 | 2.16 | 0.022 | 1.82 |
| 2 | 94 | 52.4 | 0.010 | 1.80 | 0.027 | 2.51 |
| 3 | 116 | 52.3 | 0.016 | 2.22 | 0.028 | 2.15 |
| 4 | 136 | 51.8 | 0.018 | 2.62 | 0.031 | 2.21 |
| 6 | 118 | 52.4 | 0.025 | 2.26 | 0.033 | 1.74 |
| 8 | 125 | 52.3 | 0.036 | 2.38 | 0.051 | 1.88 |
| 10 | 110 | 52.1 | 0.463 | 2.14 | 0.334 | 2.12 |
| 11 | 148 | 52.3 | 1.748 | 2.92 | 0.940 | 1.92 |
| 12 | 92 | 52.9 | 6.217 | 1.98 | 2.487 | 2.35 |
| 13 | 120 | 51.9 | 7.450 | 2.69 | 5.480 | 3.13 |
| 14 | 120 | 53.1 | 23.573 | 4.08 | 9.760 | 2.18 |
| 15 | 141 | 53.4 | 27.672 | 5.49 | 16.964 | 3.48 |
| 17 | 114 | 52.8 | 53.269 | inf | 38.339 | 6.40 |
| 20 | 89 | 52.5 | 53.330 | inf | 52.788 | inf |
| 30 | 120 | 52.2 | 52.653 | inf | 53.148 | inf |

## Cells updated per second, E. coli 536 (m = 100)

Sellers and Myers: m*n DP cells represented per pattern (Myers updates 64 per word operation); banded: cells actually computed.

| k | Sellers cells/s | banded cells computed | banded cells/s | Myers cells represented | Myers cells/s | Myers/banded cells/s | banded ms | Myers ms |
|---|---|---|---|---|---|---|---|---|
| 1 | 7.53e+08 | 3.09e+08 | 2.83e+08 | 9.78e+09 | 8.72e+09 | 31x | 54.6 | 56.0 |
| 2 | 7.78e+08 | 4.91e+08 | 3.66e+08 | 9.88e+09 | 9.50e+09 | 26x | 67.1 | 52.0 |
| 3 | 7.87e+08 | 6.79e+08 | 4.26e+08 | 9.88e+09 | 9.36e+09 | 22x | 79.7 | 52.8 |
| 4 | 8.19e+08 | 8.64e+08 | 4.52e+08 | 9.78e+09 | 9.15e+09 | 20x | 95.6 | 53.4 |
| 6 | 8.13e+08 | 1.23e+09 | 5.34e+08 | 9.88e+09 | 8.78e+09 | 16x | 115.6 | 56.3 |
| 8 | 7.71e+08 | 1.66e+09 | 5.76e+08 | 9.68e+09 | 8.83e+09 | 15x | 144.2 | 54.8 |
| 10 | 7.49e+08 | 2.00e+09 | 5.99e+08 | 9.68e+09 | 9.21e+09 | 15x | 167.4 | 52.5 |
| 11 | 8.24e+08 | 2.21e+09 | 6.31e+08 | 9.58e+09 | 8.93e+09 | 14x | 175.0 | 53.6 |
| 12 | 6.97e+08 | 2.40e+09 | 6.33e+08 | 9.48e+09 | 8.68e+09 | 14x | 189.7 | 54.6 |
| 13 | 7.85e+08 | 2.60e+09 | 6.45e+08 | 9.48e+09 | 7.91e+09 | 12x | 201.6 | 59.9 |
| 14 | 7.47e+08 | 2.78e+09 | 6.13e+08 | 9.48e+09 | 8.82e+09 | 14x | 227.0 | 53.8 |
| 15 | 7.82e+08 | 2.98e+09 | 6.38e+08 | 9.58e+09 | 8.93e+09 | 14x | 233.1 | 53.6 |
| 17 | 7.04e+08 | 3.37e+09 | 6.51e+08 | 9.38e+09 | 8.84e+09 | 14x | 259.3 | 53.1 |
| 20 | 8.16e+08 | 3.97e+09 | 6.90e+08 | 9.48e+09 | 8.93e+09 | 13x | 287.3 | 53.1 |
| 30 | 7.77e+08 | 5.83e+09 | 7.09e+08 | 9.88e+09 | 8.85e+09 | 12x | 410.7 | 55.8 |

## Stress scenarios, E. coli 536 + ~1% injected microsatellite/homopolymer arrays (m = 100, 20 reads each)

Recall by read = fraction of reads (that have a true match within k) for which the heuristic reports at least one match; recall by end = fraction of true end positions reported. Exact methods are asserted equal to Myers in every row.

| scenario | k | one-hit read | one-hit end | two-hit read | two-hit end | filter hits | filter_fa hits | filter ms | filter_fa ms | Myers ms |
|---|---|---|---|---|---|---|---|---|---|---|
| uniform | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17 | 4 | 0.026 | 0.046 | 51.9 |
| uniform | 10 | 1.00 | 1.00 | 0.95 | 0.97 | 463 | 186 | 0.796 | 0.335 | 52.0 |
| uniform | 15 | 1.00 | 1.00 | 0.45 | 1.00 | 25050 | 15693 | 24.903 | 17.881 | 53.1 |
| clustered | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 22 | 8 | 0.033 | 0.054 | 53.1 |
| clustered | 10 | 1.00 | 1.00 | 1.00 | 1.00 | 3256 | 2446 | 0.880 | 0.715 | 53.5 |
| clustered | 15 | 1.00 | 1.00 | 1.00 | 1.00 | 28859 | 18291 | 26.179 | 19.478 | 53.4 |
| indel_heavy | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17 | 5 | 0.028 | 0.047 | 53.0 |
| indel_heavy | 10 | 1.00 | 1.00 | 0.90 | 0.98 | 679 | 248 | 0.955 | 0.405 | 52.5 |
| indel_heavy | 15 | 0.95 | 0.97 | 0.60 | 0.79 | 33642 | 16342 | 31.604 | 18.600 | 51.1 |
| lowcomplex_edge | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 6878 | 12 | 0.469 | 0.046 | 50.5 |
| lowcomplex_edge | 10 | 1.00 | 1.00 | 0.95 | 0.97 | 12693 | 2839 | 1.185 | 2.415 | 51.6 |
| lowcomplex_edge | 15 | 1.00 | 1.00 | 0.65 | 0.59 | 38737 | 26362 | 22.892 | 19.427 | 52.3 |
| lowcomplex_inside | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17734 | 7111 | 0.583 | 0.665 | 46.7 |
| lowcomplex_inside | 10 | 1.00 | 1.00 | 0.95 | 0.99 | 18046 | 11129 | 1.565 | 1.183 | 50.2 |
| lowcomplex_inside | 15 | 1.00 | 0.91 | 0.75 | 0.68 | 51356 | 33393 | 18.321 | 13.643 | 49.3 |

## Results by dataset (m = 100; ms per pattern; index build once)

| dataset | length (bp) | k | Sellers | Myers | Edlib (best/read) | WFA2 (best/read) | filter | filter_fa | one-hit recall (read) | two-hit recall (read) | index build (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| E. coli 536 | 4938920 | 5 | 608.19 | 52.63 | 37.92 | 163.87 | 0.0256 | 0.0401 | 1.00 | 1.00 | 133.9 |
| E. coli 536 | 4938920 | 10 | 637.83 | 57.89 | 28.40 | 304.88 | 0.5517 | 0.3298 | 1.00 | 0.95 | 120.6 |
| SARS-CoV-2 (MT192765.1) | 29829 | 5 | 3.69 | 0.32 | 0.15 | 0.49 | 0.0041 | 0.0432 | 1.00 | 1.00 | 22.1 |
| SARS-CoV-2 (MT192765.1) | 29829 | 10 | 4.96 | 0.32 | 0.15 | 0.73 | 0.0163 | 0.0563 | 1.00 | 0.95 | 24.0 |
| phage lambda | 48502 | 5 | 6.14 | 0.51 | 0.50 | 1.14 | 0.0035 | 0.0254 | 1.00 | 1.00 | 20.3 |
| phage lambda | 48502 | 10 | 6.07 | 0.53 | 0.24 | 1.56 | 0.0153 | 0.0393 | 1.00 | 0.95 | 24.1 |
| human mtDNA | 16569 | 5 | 2.05 | 0.17 | 0.09 | 0.27 | 0.0039 | 0.0320 | 1.00 | 1.00 | 20.2 |
| human mtDNA | 16569 | 10 | 2.22 | 0.18 | 0.09 | 0.47 | 0.0082 | 0.0369 | 1.00 | 0.95 | 20.5 |

## Key numbers (main text, m = 100)

Myers, `filter`, `filter_fa`: median [10th, 90th percentile] over 100 reads, ms per query. Other columns: mean over 20 reads (Sellers 1, WFA2 5); recall columns are heuristic read recall (20 reads).

| k | k/m | Sellers | Banded | Myers | Edlib (best/read) | WFA2 (best/read) | filter | filter_fa | one-hit | two-hit | fallbacks (of 100) filter | filter_fa |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.01 | 650 | 54.6 | 51.9 [49.9, 53.7] | 28.6 | 66.4 | 0.007 [0.002, 0.017] | 0.022 [0.020, 0.026] | 0.065 (rec 1.00) | 0.034 (rec 1.00) | 0 | 0 |
| 2 | 0.02 | 635 | 67.1 | 52.4 [49.9, 57.3] | 25.4 | 87.4 | 0.010 [0.006, 0.021] | 0.027 [0.024, 0.032] | 0.056 (rec 1.00) | 0.031 (rec 1.00) | 0 | 0 |
| 3 | 0.03 | 627 | 79.7 | 52.3 [50.4, 55.1] | 24.3 | 99.8 | 0.016 [0.008, 0.027] | 0.028 [0.025, 0.033] | 0.064 (rec 1.00) | 0.031 (rec 1.00) | 0 | 0 |
| 4 | 0.04 | 597 | 95.6 | 51.8 [49.8, 54.0] | 25.5 | 121.7 | 0.018 [0.010, 0.034] | 0.031 [0.029, 0.036] | 0.054 (rec 1.00) | 0.030 (rec 1.00) | 0 | 0 |
| 6 | 0.06 | 608 | 115.6 | 52.4 [50.8, 58.0] | 25.0 | 167.7 | 0.025 [0.014, 0.045] | 0.033 [0.029, 0.040] | 0.076 (rec 1.00) | 0.030 (rec 1.00) | 0 | 0 |
| 8 | 0.08 | 628 | 144.2 | 52.3 [50.1, 56.0] | 26.3 | 215.8 | 0.036 [0.019, 0.106] | 0.051 [0.041, 0.068] | 0.051 (rec 1.00) | 0.026 (rec 0.95) | 0 | 0 |
| 10 | 0.10 | 646 | 167.4 | 52.1 [50.2, 54.0] | 31.3 | 267.1 | 0.463 [0.292, 1.631] | 0.334 [0.209, 0.479] | 0.056 (rec 1.00) | 0.022 (rec 0.95) | 0 | 0 |
| 11 | 0.11 | 581 | 175.0 | 52.3 [50.9, 53.3] | 31.2 | 252.0 | 1.748 [1.314, 2.233] | 0.940 [0.538, 1.489] | 0.040 (rec 1.00) | 0.022 (rec 0.70) | 0 | 0 |
| 12 | 0.12 | 680 | 189.7 | 52.9 [50.7, 58.5] | 28.2 | 284.2 | 6.217 [2.438, 8.089] | 2.487 [1.623, 3.874] | 0.053 (rec 1.00) | 0.026 (rec 0.75) | 0 | 0 |
| 13 | 0.13 | 604 | 201.6 | 51.9 [50.3, 53.6] | 26.7 | 331.4 | 7.450 [6.136, 21.031] | 5.480 [3.669, 8.352] | 0.047 (rec 1.00) | 0.027 (rec 0.75) | 0 | 0 |
| 14 | 0.14 | 635 | 227.0 | 53.1 [50.7, 61.7] | 28.0 | 417.4 | 23.573 [16.926, 27.406] | 9.760 [7.189, 13.653] | 0.045 (rec 1.00) | 0.025 (rec 0.70) | 0 | 0 |
| 15 | 0.15 | 613 | 233.1 | 53.4 [50.7, 59.5] | 24.9 | 318.8 | 27.672 [23.784, 36.925] | 16.964 [12.906, 22.990] | 0.047 (rec 1.00) | 0.029 (rec 0.60) | 8 | 0 |
| 17 | 0.17 | 667 | 259.3 | 52.8 [50.6, 56.8] | 24.9 | 395.9 | 53.269 [50.938, 65.751] | 38.339 [29.633, 51.766] | 0.045 (rec 0.85) | 0.020 (rec 0.45) | 99 | 7 |
| 20 | 0.20 | 581 | 287.3 | 52.5 [50.4, 57.0] | 28.6 | 431.8 | 53.330 [51.271, 59.960] | 52.788 [51.137, 56.207] | 0.041 (rec 0.95) | 0.018 (rec 0.40) | 100 | 100 |
| 30 | 0.30 | 636 | 410.7 | 52.2 [50.1, 54.2] | 32.0 | 591.5 | 52.653 [50.549, 54.026] | 53.148 [51.473, 57.623] | 0.037 (rec 0.30) | 0.018 (rec 0.00) | 100 | 100 |

## Instructions per cell and branch behaviour: banded vs Myers

Hardware performance counters: perf installed = False, CPU PMU device present = False. In this VM neither is true, so hardware counters could not be used. The numbers below come from `valgrind --tool=cachegrind --branch-sim=yes`: instruction counts are exact executed-instruction counts; branch mispredictions come from cachegrind's simulated predictor and are indicative only. They do not show cache misses or real timing.

Text: first 150,000 bp of E. coli 536; m = 100; per pattern (mean of 5). Myers cells = m * n DP cells represented (64 per word operation); banded cells = cells actually computed.

| k | algorithm | instructions per text char | instructions per cell | cells per text char | branches per text char | simulated mispredicts per 1000 chars | simulated mispredict rate |
|---|---|---|---|---|---|---|---|
| 5 | banded | 229.8 | 21.43 | 10.7 | 16.29 | 2024.7 | 12.43% |
| 5 | myers | 98.7 | 0.99 | 100.0 | 9.01 | 621.8 | 6.90% |
| 20 | banded | 754.3 | 18.92 | 39.9 | 45.42 | 1952.4 | 4.30% |
| 20 | myers | 98.8 | 0.99 | 100.0 | 9.01 | 639.4 | 7.10% |

At k = 5: banded executes 2.3x the instructions of Myers per text character while computing 9x fewer cells; instructions per cell are 21.4 (banded) against 0.99 (Myers).

At k = 20: banded executes 7.6x the instructions of Myers per text character while computing 3x fewer cells; instructions per cell are 18.9 (banded) against 0.99 (Myers).

