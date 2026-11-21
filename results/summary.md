# Benchmark summary

Main text: **E. coli 536** (see `data/SOURCES.txt`; no K-12 genome was reachable). Search time is per pattern, single thread; index build is reported separately.

Edlib and WFA2 return the BEST match per read (minimum distance, ends at that minimum), not all matches within k, and include Python call overhead; compare them on time only.

Cost-model constants (ns), calibrated on synthetic genome_synth.txt and evaluated on E. coli 536: c_hit=107.9, c_win=170, c_bp=7.52, a_scan=8.10, c_sel=15.0


Headline experiment (E. coli 536, m = 100, 100 reads per setting): values are the **median [10th, 90th percentile] over reads** of per-read query time; other experiments use means over 20 reads (Sellers 1, WFA2 5).


## Frequency-aware vs equal split (pure filter, no fallback), E. coli 536

| k | hits equal | hits fa | hit ratio (medians) | bp verified equal | bp verified fa | time equal (ms) | time fa (ms) | fa speedup vs equal | fa speedup vs Myers |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 5 [2, 10] | 1 [1, 2] | 5.0x | 410 [104, 927] | 102 [101, 103] | 0.007 [0.002, 0.017] | 0.006 [0.005, 0.007] | 1.20x | 12338x |
| 2 | 9 [5, 16] | 1 [1, 2] | 9.0x | 736 [404, 1477] | 104 [103, 208] | 0.012 [0.006, 0.025] | 0.006 [0.006, 0.008] | 1.98x | 12022x |
| 3 | 12 [6, 18] | 2 [1, 3] | 6.0x | 1041 [526, 1739] | 107 [105, 214] | 0.018 [0.009, 0.031] | 0.010 [0.007, 0.012] | 1.88x | 7784x |
| 4 | 14 [9, 25] | 3 [1, 4] | 4.5x | 1242 [641, 2424] | 110 [107, 322] | 0.021 [0.012, 0.043] | 0.008 [0.007, 0.012] | 2.58x | 8908x |
| 6 | 18 [11, 31] | 4 [3, 7] | 4.5x | 1632 [1008, 3145] | 237 [114, 666] | 0.028 [0.017, 0.052] | 0.011 [0.008, 0.018] | 2.50x | 6431x |
| 8 | 24 [14, 68] | 15 [8, 25] | 1.6x | 2437 [1381, 6990] | 1400 [708, 2442] | 0.041 [0.023, 0.136] | 0.032 [0.019, 0.053] | 1.29x | 2309x |
| 10 | 326 [216, 1184] | 212 [115, 321] | 1.5x | 38940 [26171, 134831] | 24740 [13387, 37104] | 0.651 [0.435, 2.245] | 0.445 [0.245, 0.651] | 1.46x | 164x |
| 11 | 1260 [978, 1592] | 667 [347, 1068] | 1.9x | 150050 [116381, 187210] | 80439 [43311, 125188] | 2.600 [1.995, 3.797] | 1.299 [0.727, 2.087] | 2.00x | 56x |
| 12 | 4828 [1813, 6220] | 1656 [1071, 2639] | 2.9x | 558434 [227061, 711286] | 204320 [133726, 310830] | 8.909 [3.533, 11.152] | 3.479 [2.246, 5.151] | 2.56x | 21x |
| 13 | 5494 [4530, 19394] | 3774 [2303, 6031] | 1.5x | 654841 [544943, 1873433] | 456108 [288378, 689800] | 10.379 [8.638, 29.496] | 8.068 [4.682, 11.411] | 1.29x | 9x |
| 14 | 20796 [13438, 26218] | 7370 [5305, 11344] | 2.8x | 2041166 [1413272, 2405927] | 867459 [639550, 1244160] | 35.506 [23.930, 41.833] | 14.560 [10.804, 20.360] | 2.44x | 5x |
| 15 | 23647 [19932, 28911] | 14000 [9854, 20383] | 1.7x | 2277172 [2022391, 2639100] | 1523332 [1110920, 1984815] | 36.786 [32.127, 45.955] | 24.077 [18.255, 32.216] | 1.53x | 3x |
| 17 | 97332 [83921, 111570] | 41428 [28912, 59795] | 2.3x | 4570596 [4396164, 4663902] | 3291704 [2696102, 3903638] | 74.782 [70.472, 78.689] | 56.166 [43.194, 69.117] | 1.33x | 1x |
| 20 | 435776 [128342, 492056] | 144428 [108145, 209961] | 3.0x | 4938868 [4816463, 4938916] | 4849136 [4684707, 4914649] | 113.553 [82.646, 123.957] | 82.764 [78.171, 92.561] | 1.37x | 1x |
| 30 | 2504724 [2330961, 2638079] | 1747860 [1409690, 2349176] | 1.4x | 4938920 [4938920, 4938920] | 4938920 [4938920, 4938920] | 359.733 [334.281, 408.835] | 258.324 [220.209, 340.001] | 1.39x | 0x |

## Cost model (predicted vs measured), E. coli 536

Fit: NNLS (scipy.optimize.nnls) on synthetic genome_synth.txt; evaluated here on E. coli 536 (a different text, but NOT blind: the calibration range was widened after a first fit failed on it; see the blind-evaluation section). c_win = 170 ns (unconstrained least squares gave 170, so the sign constraint was not active). Hit counts come from the index; windows and merged bp are estimated assuming uniform hit placement.

- Crossover filter vs Myers: measured k/m = 0.169 (k = 16.9), predicted k/m = 0.170 (k = 17.0)
- Crossover filter_fa vs Myers: measured k/m = 0.190 (k = 19.0), predicted k/m = 0.193 (k = 19.3)
- filter time prediction: median factor error 1.09x, worst 1.34x over 14 values of k (measured median > 0.01 ms; prediction ratio uses medians)
- filter_fa time prediction: median factor error 1.04x, worst 1.41x over 11 values of k (measured median > 0.01 ms; prediction ratio uses medians)
- Myers scan prediction: median factor error 1.09x, worst 1.11x

| k | filter measured (ms) | predicted | filter_fa measured | predicted | Myers measured | predicted |
|---|---|---|---|---|---|---|
| 1 | 0.007 [0.002, 0.017] | 0.009 | 0.006 [0.005, 0.007] | 0.006 | 73.0 [70.3, 88.8] | 80.0 |
| 2 | 0.012 [0.006, 0.025] | 0.016 | 0.006 [0.006, 0.008] | 0.008 | 74.3 [70.2, 83.0] | 80.0 |
| 3 | 0.018 [0.009, 0.031] | 0.022 | 0.010 [0.007, 0.012] | 0.010 | 74.0 [70.5, 78.7] | 80.0 |
| 4 | 0.021 [0.012, 0.043] | 0.026 | 0.008 [0.007, 0.012] | 0.011 | 72.4 [68.9, 74.4] | 80.0 |
| 6 | 0.028 [0.017, 0.052] | 0.036 | 0.011 [0.008, 0.018] | 0.016 | 72.7 [70.2, 76.0] | 80.0 |
| 8 | 0.041 [0.023, 0.136] | 0.049 | 0.032 [0.019, 0.053] | 0.037 | 73.9 [70.4, 77.0] | 80.0 |
| 10 | 0.651 [0.435, 2.245] | 0.682 | 0.445 [0.245, 0.651] | 0.461 | 72.8 [70.2, 78.2] | 80.0 |
| 11 | 2.600 [1.995, 3.797] | 2.608 | 1.299 [0.727, 2.087] | 1.412 | 72.5 [70.3, 76.0] | 80.0 |
| 12 | 8.909 [3.533, 11.152] | 9.681 | 3.479 [2.246, 5.151] | 3.550 | 72.8 [69.5, 77.6] | 80.0 |
| 13 | 10.379 [8.638, 29.496] | 11.281 | 8.068 [4.682, 11.411] | 7.901 | 74.8 [71.2, 80.5] | 80.0 |
| 14 | 35.506 [23.930, 41.833] | 35.152 | 14.560 [10.804, 20.360] | 14.868 | 74.2 [70.5, 85.6] | 80.0 |
| 15 | 36.786 [32.127, 45.955] | 39.360 | 24.077 [18.255, 32.216] | 26.336 | 73.5 [70.1, 80.5] | 80.0 |
| 17 | 74.782 [70.472, 78.689] | 80.780 | 56.166 [43.194, 69.117] | 57.187 | 73.2 [69.8, 82.2] | 80.0 |
| 20 | 113.553 [82.646, 123.957] | 121.272 | 82.764 [78.171, 92.561] | 89.066 | 72.0 [69.4, 75.9] | 80.0 |
| 30 | 359.733 [334.281, 408.835] | 344.456 | 258.324 [220.209, 340.001] | 262.848 | 73.1 [70.1, 82.5] | 80.0 |

## Index build vs query time and break-even, E. coli 536

Q* = index_ms / (Myers median ms per query - filter median ms per query); the number of queries after which building the index plus filtering beats repeated plain Myers scans (never, if the filter is not faster). Query times are medians over 100 reads (with the cost-model fallback active); the index build time is one measurement per setting.

| k | index build (ms) | Myers (ms/query) | filter (ms/query) | Q* filter | filter_fa (ms/query) | Q* filter_fa |
|---|---|---|---|---|---|---|
| 1 | 112 | 73.0 | 0.007 | 1.53 | 0.009 | 1.85 |
| 2 | 106 | 74.3 | 0.013 | 1.42 | 0.006 | 1.37 |
| 3 | 108 | 74.0 | 0.018 | 1.46 | 0.007 | 1.50 |
| 4 | 109 | 72.4 | 0.021 | 1.51 | 0.009 | 1.64 |
| 6 | 109 | 72.7 | 0.029 | 1.50 | 0.017 | 1.83 |
| 8 | 144 | 73.9 | 0.042 | 1.95 | 0.030 | 1.41 |
| 10 | 113 | 72.8 | 0.648 | 1.56 | 0.436 | 1.60 |
| 11 | 105 | 72.5 | 2.621 | 1.50 | 1.273 | 1.48 |
| 12 | 106 | 72.8 | 8.854 | 1.66 | 3.379 | 1.57 |
| 13 | 112 | 74.8 | 10.489 | 1.74 | 7.612 | 2.26 |
| 14 | 108 | 74.2 | 32.675 | 2.60 | 13.691 | 1.82 |
| 15 | 145 | 73.5 | 37.718 | 4.04 | 24.540 | 2.20 |
| 17 | 131 | 73.2 | 72.342 | 155.74 | 54.006 | 7.42 |
| 20 | 105 | 72.0 | 72.207 | inf | 72.326 | inf |
| 30 | 121 | 73.1 | 72.641 | 247.57 | 72.634 | 208.95 |

## Cost model: which figures are blind and which are fitted

Protocol (fixed before the blind text was read): constants are fitted by NNLS on the synthetic genome X (synthetic genome_synth.txt) with k up to 25, then frozen. Prediction error = median over k of the factor max(pred/measured, measured/pred) between the median predicted and the median measured per-read time of the pure filter (fallback off), for k where the measured median exceeds 0.01 ms. The fit was not changed after seeing any evaluation text.

| text | status | reads per k | filter error (median / worst) | filter_fa error (median / worst) | Myers error (median) | crossover filter measured / predicted (k/m) | crossover filter_fa measured / predicted (k/m) |
|---|---|---|---|---|---|---|---|
| Synthetic genome X | fitted text, new reads | 30 | 1.08x / 1.39x | 1.07x / 1.17x | 1.08x | 0.175 / 0.180 | 0.187 / 0.188 |
| E. coli 536 | not fitted, range chosen after seeing it | 100 | 1.09x / 1.34x | 1.04x / 1.41x | 1.09x | 0.169 / 0.170 | 0.190 / 0.193 |
| B. fragilis | blind | 100 | 1.10x / 1.27x | 1.06x / 1.19x | 1.11x | 0.169 / 0.170 | 0.186 / 0.192 |

Fitted: the constants themselves and the synthetic-text row. Not blind: the E. coli row (the calibration range was widened to k <= 25 after a first fit restricted to k <= 15 failed on this text). Blind: the B. fragilis row only.


## Cells updated per second, E. coli 536 (m = 100)

Sellers and Myers: m*n DP cells represented per pattern (Myers updates 64 per word operation); banded: cells actually computed.

| k | Sellers cells/s | banded cells computed | banded cells/s | Myers cells represented | Myers cells/s | Myers/banded cells/s | banded ms | Myers ms |
|---|---|---|---|---|---|---|---|---|
| 1 | 5.59e+08 | 3.09e+08 | 2.06e+08 | 9.78e+09 | 6.76e+09 | 33x | 75.0 | 72.3 |
| 2 | 5.71e+08 | 4.91e+08 | 2.63e+08 | 9.88e+09 | 6.52e+09 | 25x | 93.5 | 75.7 |
| 3 | 5.75e+08 | 6.79e+08 | 3.04e+08 | 9.88e+09 | 6.59e+09 | 22x | 111.6 | 75.0 |
| 4 | 5.61e+08 | 8.64e+08 | 3.29e+08 | 9.78e+09 | 6.67e+09 | 20x | 131.5 | 73.3 |
| 6 | 5.60e+08 | 1.23e+09 | 3.75e+08 | 9.88e+09 | 6.51e+09 | 17x | 164.5 | 75.9 |
| 8 | 5.33e+08 | 1.66e+09 | 4.02e+08 | 9.68e+09 | 6.35e+09 | 16x | 206.5 | 76.2 |
| 10 | 5.28e+08 | 2.00e+09 | 4.11e+08 | 9.68e+09 | 6.14e+09 | 15x | 244.1 | 78.9 |
| 11 | 5.72e+08 | 2.21e+09 | 4.37e+08 | 9.58e+09 | 6.58e+09 | 15x | 252.9 | 72.8 |
| 12 | 5.21e+08 | 2.40e+09 | 4.26e+08 | 9.48e+09 | 6.38e+09 | 15x | 282.1 | 74.3 |
| 13 | 5.71e+08 | 2.60e+09 | 4.49e+08 | 9.48e+09 | 6.38e+09 | 14x | 289.9 | 74.3 |
| 14 | 5.65e+08 | 2.78e+09 | 4.44e+08 | 9.48e+09 | 6.53e+09 | 15x | 313.2 | 72.6 |
| 15 | 5.36e+08 | 2.98e+09 | 4.62e+08 | 9.58e+09 | 6.37e+09 | 14x | 321.8 | 75.2 |
| 17 | 5.15e+08 | 3.37e+09 | 4.69e+08 | 9.38e+09 | 5.86e+09 | 12x | 359.7 | 80.0 |
| 20 | 5.53e+08 | 3.97e+09 | 4.76e+08 | 9.48e+09 | 6.55e+09 | 14x | 416.2 | 72.4 |
| 30 | 5.46e+08 | 5.83e+09 | 5.02e+08 | 9.88e+09 | 6.13e+09 | 12x | 580.8 | 80.6 |

## Stress scenarios, E. coli 536 + ~1% injected microsatellite/homopolymer arrays (m = 100, 20 reads each)

Recall by read = fraction of reads (that have a true match within k) for which the heuristic reports at least one match; recall by end = fraction of true end positions reported. Exact methods are asserted equal to Myers in every row.

| scenario | k | one-hit read | one-hit end | two-hit read | two-hit end | filter hits | filter_fa hits | filter ms | filter_fa ms | Myers ms |
|---|---|---|---|---|---|---|---|---|---|---|
| uniform | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17 | 4 | 0.033 | 0.015 | 78.0 |
| uniform | 10 | 1.00 | 1.00 | 0.95 | 0.97 | 463 | 186 | 1.192 | 0.375 | 76.8 |
| uniform | 15 | 1.00 | 1.00 | 0.45 | 1.00 | 25050 | 15693 | 34.246 | 24.592 | 79.3 |
| clustered | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 22 | 8 | 0.035 | 0.016 | 76.0 |
| clustered | 10 | 1.00 | 1.00 | 1.00 | 1.00 | 3256 | 2617 | 1.291 | 0.627 | 75.3 |
| clustered | 15 | 1.00 | 1.00 | 1.00 | 1.00 | 28859 | 18291 | 35.974 | 29.718 | 76.5 |
| indel_heavy | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17 | 5 | 0.030 | 0.013 | 72.1 |
| indel_heavy | 10 | 1.00 | 1.00 | 0.90 | 0.98 | 679 | 248 | 1.362 | 0.484 | 78.3 |
| indel_heavy | 15 | 0.95 | 0.97 | 0.60 | 0.79 | 33642 | 16342 | 42.712 | 25.989 | 72.6 |
| lowcomplex_edge | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 6878 | 12 | 0.701 | 0.038 | 68.6 |
| lowcomplex_edge | 10 | 1.00 | 1.00 | 0.95 | 0.97 | 12693 | 2992 | 2.045 | 3.648 | 71.9 |
| lowcomplex_edge | 15 | 1.00 | 1.00 | 0.65 | 0.59 | 38737 | 26362 | 33.760 | 29.325 | 69.6 |
| lowcomplex_inside | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17734 | 8666 | 1.230 | 0.626 | 65.5 |
| lowcomplex_inside | 10 | 1.00 | 1.00 | 0.95 | 0.99 | 18046 | 11163 | 2.112 | 1.462 | 65.5 |
| lowcomplex_inside | 15 | 1.00 | 0.91 | 0.75 | 0.68 | 51356 | 33393 | 25.139 | 19.000 | 77.3 |

## Results by dataset (m = 100; ms per pattern; index build once)

| dataset | length (bp) | k | Sellers | Myers | Edlib (best/read) | WFA2 (best/read) | filter | filter_fa | one-hit recall (read) | two-hit recall (read) | index build (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| E. coli 536 | 4938920 | 5 | 871.19 | 75.26 | 37.87 | 178.03 | 0.0281 | 0.0123 | 1.00 | 1.00 | 111.1 |
| E. coli 536 | 4938920 | 10 | 885.42 | 77.58 | 50.01 | 383.38 | 0.7636 | 0.5043 | 1.00 | 0.95 | 112.7 |
| SARS-CoV-2 (MT192765.1) | 29829 | 5 | 5.50 | 0.49 | 0.23 | 0.75 | 0.0045 | 0.0135 | 1.00 | 1.00 | 25.7 |
| SARS-CoV-2 (MT192765.1) | 29829 | 10 | 7.03 | 0.63 | 0.23 | 1.77 | 0.0172 | 0.0249 | 1.00 | 0.95 | 33.7 |
| phage lambda | 48502 | 5 | 12.90 | 0.72 | 0.37 | 1.51 | 0.0039 | 0.0100 | 1.00 | 1.00 | 26.6 |
| phage lambda | 48502 | 10 | 8.71 | 0.72 | 0.36 | 1.97 | 0.0191 | 0.0203 | 1.00 | 0.95 | 31.7 |
| human mtDNA | 16569 | 5 | 3.48 | 0.31 | 0.13 | 0.43 | 0.0061 | 0.0099 | 1.00 | 1.00 | 31.2 |
| human mtDNA | 16569 | 10 | 2.97 | 0.25 | 0.24 | 1.02 | 0.0098 | 0.0163 | 1.00 | 0.95 | 24.8 |

## Key numbers (main text, m = 100)

Myers, `filter`, `filter_fa`: median [10th, 90th percentile] over 100 reads, ms per query. Other columns: mean over 20 reads (Sellers 1, WFA2 5); recall columns are heuristic read recall (20 reads).

| k | k/m | Sellers | Banded | Myers | Edlib (best/read) | WFA2 (best/read) | filter | filter_fa | one-hit | two-hit | fallbacks (of 100) filter | filter_fa |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.01 | 875 | 75.0 | 73.0 [70.3, 88.8] | 43.9 | 85.0 | 0.007 [0.002, 0.017] | 0.009 [0.008, 0.010] | 0.057 (rec 1.00) | 0.040 (rec 1.00) | 0 | 0 |
| 2 | 0.02 | 865 | 93.5 | 74.3 [70.2, 83.0] | 36.4 | 87.7 | 0.013 [0.006, 0.027] | 0.006 [0.006, 0.008] | 0.053 (rec 1.00) | 0.025 (rec 1.00) | 0 | 0 |
| 3 | 0.03 | 860 | 111.6 | 74.0 [70.5, 78.7] | 36.9 | 122.5 | 0.018 [0.009, 0.031] | 0.007 [0.007, 0.009] | 0.056 (rec 1.00) | 0.030 (rec 1.00) | 0 | 0 |
| 4 | 0.04 | 871 | 131.5 | 72.4 [68.9, 74.4] | 37.1 | 173.7 | 0.021 [0.012, 0.043] | 0.009 [0.007, 0.012] | 0.049 (rec 1.00) | 0.024 (rec 1.00) | 0 | 0 |
| 6 | 0.06 | 882 | 164.5 | 72.7 [70.2, 76.0] | 45.3 | 175.7 | 0.029 [0.018, 0.053] | 0.017 [0.013, 0.025] | 0.051 (rec 1.00) | 0.023 (rec 1.00) | 0 | 0 |
| 8 | 0.08 | 908 | 206.5 | 73.9 [70.4, 77.0] | 38.1 | 260.3 | 0.042 [0.023, 0.127] | 0.030 [0.018, 0.055] | 0.050 (rec 1.00) | 0.025 (rec 0.95) | 0 | 0 |
| 10 | 0.10 | 917 | 244.1 | 72.8 [70.2, 78.2] | 38.1 | 290.9 | 0.648 [0.435, 2.255] | 0.436 [0.252, 0.736] | 0.049 (rec 1.00) | 0.030 (rec 0.95) | 0 | 0 |
| 11 | 0.11 | 837 | 252.9 | 72.5 [70.3, 76.0] | 37.1 | 317.9 | 2.621 [2.057, 3.496] | 1.273 [0.688, 1.990] | 0.044 (rec 1.00) | 0.021 (rec 0.70) | 0 | 0 |
| 12 | 0.12 | 910 | 282.1 | 72.8 [69.5, 77.6] | 37.1 | 325.3 | 8.854 [3.397, 11.238] | 3.379 [2.140, 5.155] | 0.050 (rec 1.00) | 0.022 (rec 0.75) | 0 | 0 |
| 13 | 0.13 | 830 | 289.9 | 74.8 [71.2, 80.5] | 37.3 | 335.6 | 10.489 [8.730, 29.967] | 7.612 [4.766, 11.452] | 0.048 (rec 1.00) | 0.020 (rec 0.75) | 0 | 0 |
| 14 | 0.14 | 839 | 313.2 | 74.2 [70.5, 85.6] | 52.4 | 409.9 | 32.675 [22.759, 38.873] | 13.691 [10.051, 19.292] | 0.045 (rec 1.00) | 0.028 (rec 0.70) | 0 | 0 |
| 15 | 0.15 | 894 | 321.8 | 73.5 [70.1, 80.5] | 41.7 | 396.1 | 37.718 [32.249, 51.190] | 24.540 [18.134, 32.815] | 0.049 (rec 1.00) | 0.025 (rec 0.60) | 2 | 0 |
| 17 | 0.17 | 910 | 359.7 | 73.2 [69.8, 82.2] | 36.6 | 475.0 | 72.342 [69.068, 75.908] | 54.006 [41.486, 64.468] | 0.043 (rec 0.85) | 0.022 (rec 0.45) | 64 | 0 |
| 20 | 0.20 | 857 | 416.2 | 72.0 [69.4, 75.9] | 44.1 | 539.3 | 72.207 [69.524, 76.246] | 72.326 [69.364, 76.221] | 0.047 (rec 0.95) | 0.022 (rec 0.40) | 100 | 97 |
| 30 | 0.30 | 904 | 580.8 | 73.1 [70.1, 82.5] | 57.1 | 730.9 | 72.641 [69.536, 81.448] | 72.634 [69.952, 75.150] | 0.043 (rec 0.30) | 0.015 (rec 0.00) | 100 | 100 |

## Instructions per cell and branch behaviour: banded vs Myers

Hardware performance counters: perf installed = False, CPU PMU device present = False. In this VM neither is true, so hardware counters could not be used. The numbers below come from `valgrind --tool=cachegrind --branch-sim=yes`: instruction counts are exact executed-instruction counts; branch mispredictions come from cachegrind's simulated predictor and are indicative only. They do not show cache misses or real timing.

Text: first 150,000 bp of E. coli 536; m = 100; per pattern (mean of 5). Myers cells = m * n DP cells represented (64 per word operation); banded cells = cells actually computed.

| k | algorithm | instructions per text char | instructions per cell | cells per text char | branches per text char | simulated mispredicts per 1000 chars | simulated mispredict rate |
|---|---|---|---|---|---|---|---|
| 5 | banded | 229.8 | 21.43 | 10.7 | 16.29 | 2024.7 | 12.43% |
| 5 | myers | 103.7 | 1.04 | 100.0 | 9.01 | 621.8 | 6.90% |
| 20 | banded | 754.3 | 18.92 | 39.9 | 45.42 | 1952.4 | 4.30% |
| 20 | myers | 103.8 | 1.04 | 100.0 | 9.01 | 639.4 | 7.10% |

At k = 5: banded executes 2.2x the instructions of Myers per text character while computing 9x fewer cells; instructions per cell are 21.4 (banded) against 1.04 (Myers).

At k = 20: banded executes 7.3x the instructions of Myers per text character while computing 3x fewer cells; instructions per cell are 18.9 (banded) against 1.04 (Myers).


## Expected-case analysis: formula versus measured seed hits (equal split, m = 100)

`hits` = seed hits per read, mean over reads. (A) random uniform patterns: formula (k+1) N 4^-s. (B) reads sampled from the text with exactly k edits: uniform-text formula (k+1)(lambda_s + pi) and spectrum-corrected formula (k+1)[pi C_s + (1-pi) N 4^-s], pi = (1 - s/m')^k, both averaged over the read-length distribution m' = m + insertions - deletions (piece-length detail).


### Uniform random text (no repeats) (n = 4,600,000)

| k | q | s | (A) formula | (A) measured | (A) ratio | (B) uniform formula | (B) spectrum-corrected | (B) measured | measured / uniform | measured / corrected | C_s | 1 + lambda_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 50 | 11 | 2.2 | 2.3 | 1.04 | 4.0 | 4.0 | 4.0 | 0.99 | 0.99 | 2.10 | 2.10 |
| 2 | 33 | 11 | 3.3 | 3.3 | 1.00 | 5.7 | 5.7 | 5.8 | 1.02 | 1.02 | 2.10 | 2.10 |
| 3 | 25 | 11 | 4.4 | 4.2 | 0.97 | 7.2 | 7.2 | 7.6 | 1.06 | 1.06 | 2.10 | 2.10 |
| 4 | 20 | 11 | 5.5 | 5.5 | 1.00 | 8.6 | 8.6 | 8.6 | 0.99 | 0.99 | 2.10 | 2.10 |
| 6 | 14 | 11 | 7.7 | 7.5 | 0.98 | 11.2 | 11.2 | 11.2 | 1.00 | 1.00 | 2.10 | 2.10 |
| 8 | 11 | 11 | 9.9 | 9.7 | 0.99 | 21.2 | 21.2 | 20.8 | 0.98 | 0.98 | 2.10 | 2.10 |
| 10 | 9 | 9 | 193.0 | 193.8 | 1.00 | 361.2 | 361.2 | 350.6 | 0.97 | 0.97 | 18.54 | 18.55 |
| 11 | 8 | 8 | 842.3 | 840.6 | 1.00 | 966.1 | 966.0 | 960.0 | 0.99 | 0.99 | 71.18 | 71.19 |
| 12 | 7 | 7 | 3649.9 | 3649.5 | 1.00 | 3359.8 | 3359.7 | 3331.3 | 0.99 | 0.99 | 281.73 | 281.76 |
| 13 | 7 | 7 | 3930.7 | 3921.1 | 1.00 | 6287.8 | 6287.6 | 5879.0 | 0.93 | 0.94 | 281.73 | 281.76 |
| 14 | 6 | 6 | 16845.7 | 16839.5 | 1.00 | 15971.2 | 15970.9 | 16444.7 | 1.03 | 1.03 | 1123.98 | 1124.05 |
| 15 | 6 | 6 | 17968.7 | 17929.1 | 1.00 | 22145.6 | 22145.1 | 22475.7 | 1.01 | 1.01 | 1123.98 | 1124.05 |
| 17 | 5 | 5 | 80859.3 | 80802.0 | 1.00 | 80286.4 | 80285.2 | 80881.5 | 1.01 | 1.01 | 4493.02 | 4493.18 |
| 20 | 4 | 4 | 377343.5 | 377403.1 | 1.00 | 346425.7 | 346424.1 | 344398.7 | 0.99 | 0.99 | 17969.56 | 17969.74 |

Predicted crossover with a Myers scan (cost-model constants, formula hits): uniform-text formula k = 17.4 (k/m = 0.174); spectrum-corrected formula k = 17.4 (k/m = 0.174).

### Synthetic genome (3% injected repeats) (n = 4,600,000)

| k | q | s | (A) formula | (A) measured | (A) ratio | (B) uniform formula | (B) spectrum-corrected | (B) measured | measured / uniform | measured / corrected | C_s | 1 + lambda_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 50 | 11 | 2.2 | 2.1 | 0.97 | 4.0 | 4.1 | 4.1 | 1.02 | 1.00 | 2.14 | 2.10 |
| 2 | 33 | 11 | 3.3 | 3.3 | 1.01 | 5.7 | 5.8 | 5.9 | 1.04 | 1.03 | 2.14 | 2.10 |
| 3 | 25 | 11 | 4.4 | 4.3 | 0.97 | 7.2 | 7.3 | 7.2 | 1.00 | 0.98 | 2.14 | 2.10 |
| 4 | 20 | 11 | 5.5 | 5.5 | 1.00 | 8.6 | 8.8 | 8.6 | 1.00 | 0.99 | 2.14 | 2.10 |
| 6 | 14 | 11 | 7.7 | 7.8 | 1.02 | 11.2 | 11.3 | 11.6 | 1.04 | 1.02 | 2.14 | 2.10 |
| 8 | 11 | 11 | 9.9 | 10.3 | 1.04 | 21.2 | 21.4 | 21.0 | 0.99 | 0.98 | 2.14 | 2.10 |
| 10 | 9 | 9 | 193.0 | 192.6 | 1.00 | 361.2 | 361.4 | 349.5 | 0.97 | 0.97 | 18.59 | 18.55 |
| 11 | 8 | 8 | 842.3 | 844.4 | 1.00 | 966.1 | 966.3 | 964.4 | 1.00 | 1.00 | 71.24 | 71.19 |
| 12 | 7 | 7 | 3649.9 | 3649.9 | 1.00 | 3359.8 | 3360.0 | 3325.1 | 0.99 | 0.99 | 281.80 | 281.76 |
| 13 | 7 | 7 | 3930.7 | 3922.9 | 1.00 | 6287.8 | 6288.0 | 5907.3 | 0.94 | 0.94 | 281.80 | 281.76 |
| 14 | 6 | 6 | 16845.7 | 16863.2 | 1.00 | 15971.2 | 15971.4 | 16400.4 | 1.03 | 1.03 | 1124.08 | 1124.05 |
| 15 | 6 | 6 | 17968.7 | 17996.8 | 1.00 | 22145.6 | 22145.8 | 22492.0 | 1.02 | 1.02 | 1124.08 | 1124.05 |
| 17 | 5 | 5 | 80859.3 | 80824.0 | 1.00 | 80286.4 | 80286.4 | 80822.7 | 1.01 | 1.01 | 4493.18 | 4493.18 |
| 20 | 4 | 4 | 377343.5 | 377303.2 | 1.00 | 346425.7 | 346425.5 | 344371.9 | 0.99 | 0.99 | 17969.72 | 17969.74 |

Predicted crossover with a Myers scan (cost-model constants, formula hits): uniform-text formula k = 17.4 (k/m = 0.174); spectrum-corrected formula k = 17.4 (k/m = 0.174).

### E. coli 536 (real) (n = 4,938,920)

| k | q | s | (A) formula | (A) measured | (A) ratio | (B) uniform formula | (B) spectrum-corrected | (B) measured | measured / uniform | measured / corrected | C_s | 1 + lambda_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 50 | 11 | 2.4 | 2.5 | 1.08 | 4.1 | 6.7 | 6.9 | 1.67 | 1.03 | 3.60 | 2.18 |
| 2 | 33 | 11 | 3.5 | 3.7 | 1.05 | 5.9 | 9.3 | 9.6 | 1.63 | 1.04 | 3.60 | 2.18 |
| 3 | 25 | 11 | 4.7 | 4.7 | 1.00 | 7.5 | 11.5 | 12.8 | 1.70 | 1.11 | 3.60 | 2.18 |
| 4 | 20 | 11 | 5.9 | 5.6 | 0.95 | 9.0 | 13.5 | 15.3 | 1.69 | 1.13 | 3.60 | 2.18 |
| 6 | 14 | 11 | 8.2 | 7.7 | 0.94 | 11.7 | 16.7 | 19.1 | 1.63 | 1.14 | 3.60 | 2.18 |
| 8 | 11 | 11 | 10.6 | 10.1 | 0.95 | 22.5 | 30.6 | 34.7 | 1.54 | 1.14 | 3.60 | 2.18 |
| 10 | 9 | 9 | 207.2 | 203.6 | 0.98 | 387.5 | 485.2 | 536.9 | 1.39 | 1.11 | 33.41 | 19.84 |
| 11 | 8 | 8 | 904.3 | 878.3 | 0.97 | 1036.9 | 1261.7 | 1356.1 | 1.31 | 1.07 | 118.63 | 76.36 |
| 12 | 7 | 7 | 3918.8 | 3967.6 | 1.01 | 3607.0 | 4254.5 | 4639.1 | 1.29 | 1.09 | 431.87 | 302.45 |
| 13 | 7 | 7 | 4220.3 | 4370.7 | 1.04 | 6750.7 | 7779.3 | 7931.6 | 1.17 | 1.02 | 431.87 | 302.45 |
| 14 | 6 | 6 | 18086.8 | 18492.5 | 1.02 | 17147.5 | 19422.4 | 21495.2 | 1.25 | 1.11 | 1587.30 | 1206.79 |
| 15 | 6 | 6 | 19292.6 | 19580.3 | 1.01 | 23776.8 | 26580.8 | 28771.2 | 1.21 | 1.08 | 1587.30 | 1206.79 |
| 17 | 5 | 5 | 86816.9 | 85768.4 | 0.99 | 86201.2 | 94052.5 | 96206.8 | 1.12 | 1.02 | 5876.23 | 4824.16 |
| 20 | 4 | 4 | 405145.5 | 408636.8 | 1.01 | 371949.1 | 394515.5 | 411230.2 | 1.11 | 1.04 | 21929.96 | 19293.64 |

Predicted crossover with a Myers scan (cost-model constants, formula hits): uniform-text formula k = 17.4 (k/m = 0.174); spectrum-corrected formula k = 17.0 (k/m = 0.170).

### Agreement (factor error over k = 1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 14, 15, 17, 20: median / worst)

| text | (A) random patterns vs formula | (B) reads vs uniform formula | (B) reads vs spectrum-corrected formula |
|---|---|---|---|
| Uniform random text (no repeats) | 1.00x / 1.04x | 1.01x / 1.07x | 1.01x / 1.07x |
| Synthetic genome (3% injected repeats) | 1.00x / 1.04x | 1.02x / 1.06x | 1.02x / 1.06x |
| E. coli 536 (real) | 1.03x / 1.08x | 1.39x / 1.70x | 1.09x / 1.14x |


## Real divergence: mt_orang segments against mt_human, and the reverse

Both are real homologous mitochondrial genomes (human, 16,569 bp; orangutan, 16,499 bp), so the edits are real evolutionary differences, not simulated ones. Segments are random substrings of the query genome; the text is the other genome (forward strand). The true edit distance is the minimum over all text positions, computed with Myers and verified with Sellers.


### mt_orang segments -> mt_human text, m = 100 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 15 (15.0% of m), 10th to 90th percentile 8 to 22 (8.0% to 22.0%), min 1, max 39; 0 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 39 | 0.2431 | 0.1190 | 0.9618 | 0.0071 | 0.0118 | 5 | 2 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 110 | 0.2543 | 0.1569 | 0.9378 | 0.1605 | 0.1172 | 135 | 73 | 0 / 200 | 1.00 | 0.96 | 1.00 | 0.98 |
| 20 | 177 | 0.2459 | 0.1273 | 0.9236 | 0.3601 | 0.2934 | 2082 | 680 | 200 / 200 | 0.99 | 0.85 | 1.00 | 0.93 |
| 25 | 190 | 0.2471 | 0.1307 | 0.9112 | 0.8607 | 0.5221 | 9006 | 2786 | 200 / 200 | 0.99 | 0.82 | 1.00 | 0.89 |
| 30 | 195 | 0.2456 | 0.1604 | 0.9181 | 1.0406 | 0.8045 | 10743 | 7682 | 200 / 200 | 0.98 | 0.80 | 0.99 | 0.86 |
| 40 | 200 | 0.2417 | 0.1545 | 0.9036 | 4.2709 | 2.5036 | 51031 | 28845 | 200 / 200 | 0.96 | 0.79 | 0.77 | 0.65 |

### mt_orang segments -> mt_human text, m = 200 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 29 (14.5% of m), 10th to 90th percentile 20 to 41 (10.0% to 20.5%), min 4, max 93; 0 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 3 | 0.6268 | 0.1613 | 1.8116 | 0.0074 | 0.0171 | 2 | 0 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 11 | 0.4885 | 0.1632 | 2.1795 | 0.0088 | 0.0192 | 4 | 1 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 20 | 22 | 0.4802 | 0.1523 | 1.8059 | 0.0396 | 0.0395 | 10 | 3 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 25 | 60 | 0.4858 | 0.2520 | 2.4985 | 0.2945 | 0.0973 | 72 | 16 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 30 | 129 | 0.4853 | 0.1681 | 1.7462 | 0.4880 | 0.4264 | 265 | 105 | 12 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 40 | 179 | 0.5078 | 0.2914 | 2.2528 | 0.7304 | 0.5902 | 4073 | 1133 | 200 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |

### mt_human segments -> mt_orang text, m = 100 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 15 (15.0% of m), 10th to 90th percentile 8 to 24 (8.0% to 24.0%), min 0, max 40; 1 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 43 | 0.2441 | 0.2104 | 1.5276 | 0.0084 | 0.0173 | 5 | 2 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 111 | 0.2411 | 0.1227 | 0.9627 | 0.1592 | 0.1397 | 136 | 69 | 0 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |
| 20 | 163 | 0.2492 | 0.1334 | 1.0017 | 0.3659 | 0.2879 | 2067 | 664 | 200 / 200 | 1.00 | 0.86 | 1.00 | 0.94 |
| 25 | 182 | 0.2398 | 0.2286 | 1.0104 | 0.9091 | 0.4236 | 8979 | 2746 | 200 / 200 | 0.99 | 0.79 | 1.00 | 0.89 |
| 30 | 188 | 0.3034 | 0.1945 | 1.0562 | 0.9918 | 0.7728 | 10736 | 7541 | 200 / 200 | 0.99 | 0.78 | 1.00 | 0.85 |
| 40 | 200 | 0.2480 | 0.1738 | 1.0041 | 4.4636 | 2.4549 | 50822 | 28455 | 200 / 200 | 0.95 | 0.73 | 0.77 | 0.63 |

### mt_human segments -> mt_orang text, m = 200 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 29 (14.5% of m), 10th to 90th percentile 21 to 41 (10.5% to 20.5%), min 6, max 85; 0 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 3 | 0.4827 | 0.2092 | 2.0635 | 0.0110 | 0.0143 | 2 | 0 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 11 | 0.5992 | 0.1195 | 2.0005 | 0.0092 | 0.0168 | 4 | 1 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 20 | 19 | 0.6040 | 0.1188 | 1.7867 | 0.0376 | 0.0312 | 10 | 3 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 25 | 47 | 0.6259 | 0.1394 | 2.1006 | 0.2908 | 0.0959 | 72 | 15 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 30 | 115 | 0.4822 | 0.2896 | 1.8552 | 0.5112 | 0.4087 | 264 | 99 | 6 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 40 | 178 | 0.4919 | 0.1743 | 2.2775 | 0.7323 | 0.5867 | 4079 | 1099 | 200 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |

Times are mean ms per segment on a 16.5 kbp text (all exact methods agree with Sellers on every segment and k; Edlib and WFA2 return the best match per read and include Python overhead). With a text this small the index costs about 20 ms to build (clearing the 4^11 table), so the filter is only worthwhile over many reads. "pure" = fallback disabled; the fallback column counts how many segments the cost-model rule sent to a full Myers scan.


## Alignment output validation

3200 random reads in 320 runs (texts of 200 to 20,000 bp, some with N; k in 0..12; uniform, indel-heavy and clustered edits; 2298 mapped, 902 unmapped): 0 failures. Checks per mapped read: distance equals the Sellers minimum (Sellers at k = d matches, at k = d - 1 does not); CIGAR edit count (X + I + D) equals the reported distance; the CIGAR replays correctly against read and text; the end is a Sellers end position. Unmapped reads have Sellers minimum distance above k.



### Worked example: `scripts/map_reads` on E. coli 536

A 60 bp stretch of E. coli 536 (positions 2000000..2000060) was copied and given 2 substitutions, 1 insertion and 1 deletion (4 edits), then mapped with k = 8:

```
scripts/map_reads data/raw/ecoli536.fa reads.fa -k 8
```
TSV output:
```
read	length	status	distance	start	end	strand	cigar	best_ends
example_read	60	mapped	4	2000000	2000060	+	10=1X15=1I14=1D12=1X6=	10
unmappable_read	60	unmapped	*	*	*	*	*	0
```
Alignment reconstructed from the CIGAR (reference on top, read below; `x` mismatch, `-` gap):
```
ATATGGCAAAAGCGCTCAGGGCGGGA-TCATCAACATCGTCACCCAGCAGCCGGACAGCAC
||||||||||x||||||||||||||| |||||||||||||| ||||||||||||x||||||
ATATGGCAAACGCGCTCAGGGCGGGAGTCATCAACATCGTC-CCCAGCAGCCGGTCAGCAC
```
The reported distance is 4, equal to the number of X, I and D operations in `10=1X15=1I14=1D12=1X6=`, and the start (2000000) matches the planted position (2000000). The second read has no alignment within k = 8 and is reported as unmapped. SAM-like output (`--format sam`):
```
@HD	VN:1.6	SO:unsorted
@SQ	SN:gi|110640213|ref|NC_008253.1|	LN:4938920
@PG	ID:map_reads	PN:map_reads
example_read	0	gi|110640213|ref|NC_008253.1|	2000001	255	10=1X15=1I14=1D12=1X6=	*	0	0	ATATGGCAAACGCGCTCAGGGCGGGAGTCATCAACATCGTCCCCAGCAGCCGGTCAGCAC	*	NM:i:4
unmappable_read	4	*	0	0	*	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	*
```


## Seed-selection objective: hits versus hits + lambda * windows

Pure filters (fallback off), mean over reads; windows = merged windows verified; lambda = (c_win + c_bp B (m + 2k)) / c_hit with the fitted constants (c_hit = 107.9, c_win = 170, c_bp = 7.52 ns), e.g. lambda = 18.3 at m = 100, k = 10. All methods return the same matches as Myers (asserted).

### Stress scenarios (E. coli 536 + about 1% microsatellite/homopolymer arrays, m = 100, 20 reads)

| scenario | k | algorithm | hits | windows | bp verified | query ms | seed selection ms |
|---|---|---|---|---|---|---|---|
| uniform | 6 | equal split | 17 | 14 | 1590 | 0.030 | 0.000 |
| uniform | 6 | fa: min hits | 4 | 2 | 241 | 0.013 | 0.007 |
| uniform | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 4 | 2 | 241 | 0.061 | 0.055 |
| uniform | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 4 | 2 | 241 | 0.072 | 0.064 |
| uniform | 10 | equal split | 463 | 448 | 54540 | 1.134 | 0.000 |
| uniform | 10 | fa: min hits | 186 | 181 | 22008 | 0.497 | 0.019 |
| uniform | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 186 | 181 | 22008 | 0.554 | 0.076 |
| uniform | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 186 | 181 | 22008 | 0.433 | 0.065 |
| uniform | 15 | equal split | 25050 | 11833 | 2165625 | 34.073 | 0.000 |
| uniform | 15 | fa: min hits | 15693 | 9462 | 1522426 | 24.105 | 0.026 |
| uniform | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16012 | 9568 | 1548422 | 25.476 | 0.154 |
| uniform | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 15693 | 9462 | 1522426 | 24.163 | 0.114 |
| clustered | 6 | equal split | 22 | 16 | 1793 | 0.034 | 0.000 |
| clustered | 6 | fa: min hits | 8 | 3 | 384 | 0.019 | 0.009 |
| clustered | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 8 | 3 | 384 | 0.059 | 0.051 |
| clustered | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 8 | 3 | 384 | 0.083 | 0.071 |
| clustered | 10 | equal split | 3256 | 513 | 62756 | 1.528 | 0.000 |
| clustered | 10 | fa: min hits | 2617 | 229 | 28295 | 0.594 | 0.015 |
| clustered | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 2616 | 228 | 28192 | 0.761 | 0.189 |
| clustered | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 2446 | 376 | 47397 | 1.560 | 0.082 |
| clustered | 15 | equal split | 28859 | 11635 | 2235284 | 35.344 | 0.000 |
| clustered | 15 | fa: min hits | 18291 | 9384 | 1530393 | 27.678 | 0.032 |
| clustered | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 19064 | 9771 | 1618158 | 26.104 | 0.222 |
| clustered | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 18291 | 9384 | 1530393 | 25.551 | 0.116 |
| indel_heavy | 6 | equal split | 17 | 14 | 1528 | 0.029 | 0.000 |
| indel_heavy | 6 | fa: min hits | 5 | 2 | 275 | 0.015 | 0.008 |
| indel_heavy | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 5 | 2 | 264 | 0.055 | 0.047 |
| indel_heavy | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 5 | 2 | 264 | 0.066 | 0.059 |
| indel_heavy | 10 | equal split | 679 | 655 | 78334 | 1.704 | 0.000 |
| indel_heavy | 10 | fa: min hits | 248 | 242 | 28661 | 0.486 | 0.013 |
| indel_heavy | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 248 | 242 | 28661 | 0.541 | 0.062 |
| indel_heavy | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 248 | 242 | 28661 | 0.532 | 0.061 |
| indel_heavy | 15 | equal split | 33642 | 11860 | 2568773 | 41.811 | 0.000 |
| indel_heavy | 15 | fa: min hits | 16342 | 10085 | 1668047 | 28.375 | 0.028 |
| indel_heavy | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16342 | 10085 | 1668047 | 27.986 | 0.125 |
| indel_heavy | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 16342 | 10085 | 1668047 | 26.982 | 0.119 |
| lowcomplex_edge | 6 | equal split | 6878 | 21 | 12819 | 0.525 | 0.000 |
| lowcomplex_edge | 6 | fa: min hits | 12 | 10 | 1148 | 0.031 | 0.008 |
| lowcomplex_edge | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 12 | 10 | 1148 | 0.699 | 0.678 |
| lowcomplex_edge | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 12 | 10 | 1148 | 0.073 | 0.053 |
| lowcomplex_edge | 10 | equal split | 12693 | 531 | 75968 | 1.619 | 0.000 |
| lowcomplex_edge | 10 | fa: min hits | 2992 | 1411 | 178281 | 3.568 | 0.017 |
| lowcomplex_edge | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 5662 | 295 | 42867 | 1.360 | 0.541 |
| lowcomplex_edge | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 2859 | 1515 | 191528 | 3.250 | 0.069 |
| lowcomplex_edge | 15 | equal split | 38737 | 10889 | 1938189 | 38.332 | 0.000 |
| lowcomplex_edge | 15 | fa: min hits | 26362 | 9964 | 1696660 | 27.670 | 0.024 |
| lowcomplex_edge | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 28203 | 10940 | 1934724 | 34.983 | 0.444 |
| lowcomplex_edge | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 26362 | 9964 | 1696660 | 27.883 | 0.116 |
| lowcomplex_inside | 6 | equal split | 17734 | 14 | 11878 | 0.830 | 0.000 |
| lowcomplex_inside | 6 | fa: min hits | 8666 | 13 | 11333 | 0.912 | 0.009 |
| lowcomplex_inside | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 7721 | 35 | 13891 | 3.121 | 2.298 |
| lowcomplex_inside | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 7133 | 121 | 23851 | 1.171 | 0.127 |
| lowcomplex_inside | 10 | equal split | 18046 | 535 | 85941 | 2.086 | 0.000 |
| lowcomplex_inside | 10 | fa: min hits | 11163 | 436 | 64069 | 1.451 | 0.012 |
| lowcomplex_inside | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 11922 | 334 | 51908 | 2.026 | 0.748 |
| lowcomplex_inside | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 11150 | 497 | 71536 | 1.755 | 0.109 |
| lowcomplex_inside | 15 | equal split | 51356 | 7508 | 1461597 | 26.797 | 0.000 |
| lowcomplex_inside | 15 | fa: min hits | 33393 | 6438 | 1121207 | 19.244 | 0.020 |
| lowcomplex_inside | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 36003 | 7417 | 1349036 | 31.724 | 0.833 |
| lowcomplex_inside | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 33415 | 6472 | 1125427 | 19.214 | 0.105 |

### Unmodified E. coli 536, uniform edits, m = 100, 100 reads (regression check)

| k | algorithm | hits | windows | bp verified | query ms | seed selection ms |
|---|---|---|---|---|---|---|
| 4 | equal split | 16 | 13 | 1478 | 0.026 | 0.000 |
| 4 | fa: min hits | 3 | 2 | 177 | 0.010 | 0.005 |
| 4 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 3 | 2 | 177 | 0.047 | 0.043 |
| 4 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 3 | 2 | 177 | 0.045 | 0.042 |
| 8 | equal split | 34 | 30 | 3479 | 0.064 | 0.000 |
| 8 | fa: min hits | 16 | 13 | 1515 | 0.034 | 0.007 |
| 8 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16 | 13 | 1509 | 0.079 | 0.051 |
| 8 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 16 | 13 | 1509 | 0.076 | 0.048 |
| 10 | equal split | 493 | 479 | 57545 | 0.944 | 0.000 |
| 10 | fa: min hits | 215 | 210 | 25147 | 0.437 | 0.013 |
| 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 215 | 210 | 25135 | 0.528 | 0.071 |
| 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 215 | 210 | 25135 | 0.487 | 0.065 |
| 12 | equal split | 4654 | 4050 | 538390 | 8.482 | 0.000 |
| 12 | fa: min hits | 1760 | 1665 | 211395 | 3.381 | 0.016 |
| 12 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 1760 | 1665 | 211377 | 3.509 | 0.090 |
| 12 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 1760 | 1665 | 211377 | 3.449 | 0.086 |
| 15 | equal split | 28807 | 12005 | 2422582 | 38.671 | 0.000 |
| 15 | fa: min hits | 14685 | 9594 | 1551347 | 25.591 | 0.027 |
| 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 14685 | 9594 | 1551347 | 24.799 | 0.120 |
| 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 14685 | 9594 | 1551347 | 24.942 | 0.121 |

**Before/after on the case that verified more bases (lowcomplex_edge, k = 10):** equal split 12693 hits, 531 windows, 75968 bp, 1.62 ms; min-hits DP 2992 hits, 1411 windows, 178281 bp, 3.57 ms; hits + lambda * windows DP 5662 hits, 295 windows, 42867 bp, 1.36 ms (seed selection 0.54 ms); same with hit-list cap 1000: 2859 hits, 1515 windows, 191528 bp, 3.25 ms (selection 0.07 ms).
`filter_fa_w` versus `filter_fa` over all 20 settings: more than 5% faster in 2, more than 5% slower in 15, within 5% in 3.
`filter_fa_w1k` versus `filter_fa` over all 20 settings: more than 5% faster in 3, more than 5% slower in 11, within 5% in 6.


## Peak memory per method (measured)

E. coli 536 (n = 4,938,920), m = 100, k = 10, 5 reads; peak resident set size (VmHWM). working memory = peak RSS minus the steady RSS after the text and reads are loaded and the heap trimmed (8.5 MB); methods whose peak does not exceed the transient peak of loading the file (about 20 MB, the file is read line by line into a string) are shown as negligible. Edlib and WFA2 run in a Python process, so their baseline differs.

| method | peak RSS (MB) | working memory of the method (MB) | expected |
|---|---|---|---|
| noop | 20.2 | negligible (below the loading peak) | text only (n bytes) plus read buffers |
| sellers | 20.3 | negligible (below the loading peak) | O(m) |
| banded | 20.3 | negligible (below the loading peak) | O(m) |
| myers | 20.3 | negligible (below the loading peak) | O(B sigma) |
| myers_pruned | 20.3 | negligible (below the loading peak) | O(B sigma) |
| lv | 20.3 | negligible (below the loading peak) | packed text n/4 bytes (1.2 MB) + L1-sized chunk buffers |
| lv_direct | 50.9 | 42.4 | 2 layers of n + 2k ints plus n flags (about 44 MB) |
| filter | 69.4 | 60.9 | index 4(4^11+1) + 4n bytes (37 MB) plus build temporaries (counter copy 17 MB, flags 5 MB) |
| filter_auto | 69.3 | 60.8 | same index |
| filter_fa | 69.3 | 60.9 | same index |
| filter_fa_w | 69.4 | 60.9 | same index |
| blast1 | 69.3 | 60.8 | same index |
| blast2 | 69.4 | 60.9 | same index |
| map | 69.4 | 60.9 | same index + O(m (m + k)) per cluster |
| edlib (Python process) | 39.3 | negligible (below the loading peak) | library; measured as ru_maxrss increase over loading the text |


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


## Phase 3A, item 2: faster Landau-Vishkin

E. coli 536, m = 100, 5 reads, best of 3 runs, ms per read. Every variant returns the same matches as Sellers (tests/test_lv.py: 600 comparisons on tandem-repeat, low-complexity and N texts; crosscheck: 300 random cases; real data: 864 comparisons; the scalar path is also tested with APPROX_LVSIMD=0).

| k | phase-2 LV | lv_direct (control) | lv_stream | lv scalar blocked | lv AVX-512 (default) | speedup of default vs phase-2 | Myers |
|---|---|---|---|---|---|---|---|
| 1 | 94.1 | 94.6 | 66.7 | 46.5 | 16.5 | 5.70x | 70.2 |
| 5 | 239.6 | 240.0 | 165.8 | 135.2 | 43.6 | 5.49x | 70.0 |
| 10 | 444.0 | 443.9 | 280.6 | 244.6 | 79.9 | 5.56x | 73.1 |
| 30 | 1232.2 | 1245.9 | 758.2 | 695.0 | 221.5 | 5.56x | 72.5 |

Variants and what they showed: the packed XOR + ctz extension alone (lv_stream) gave about 1.5x; temporal blocking (a chunk of diagonals carried through all layers in L1) gave a small and inconsistent gain (the scalar cell is instruction-bound, about 30 instructions or 4.4 ns, not memory-bound); the AVX-512 cell loop (8 cells per iteration, two gathers, lzcnt) gave the large step. A CPU without AVX-512 uses the scalar blocked loop (about 1.7x to 2.3x faster than phase 2). Texts or patterns containing N use the scalar loop.


## Phase 3A, item 3: single-word fast path for Myers and Myers with pruning

E. coli 536, 10 reads, best of 5 runs, ms per read. 'phase 2' = binary built from the final phase-2 commit; 'now' = current binary. Rows with m <= 60 are one block for every read (edits change a read's length by at most k, so every read stays within 64 bases); an earlier version of this table used m = 64, where about half the reads grow past 64 bases and need a second block, which mixed one- and two-block scans. Exactness: crosscheck (300 cases), real data (864 comparisons), tests/test_pruned.py (320 comparisons, m up to 600, k up to 200): 0 failures.

| m (blocks) | k | Myers phase 2 | Myers now | pruned phase 2 | pruned now | pruned now / Myers now | pruned now vs phase-2 pruned | Edlib (best/read) |
|---|---|---|---|---|---|---|---|---|
| 20 (1) | 2 | 41.2 | 20.6 | 50.2 | 20.6 | 1.00 | 2.43x | 36.5 |
| 32 (1) | 5 | 41.1 | 20.6 | 50.9 | 20.6 | 1.00 | 2.47x | 38.2 |
| 50 (1) | 2 | 41.1 | 20.7 | 50.1 | 20.6 | 0.99 | 2.43x | 37.7 |
| 50 (1) | 5 | 41.7 | 20.6 | 51.7 | 20.6 | 1.00 | 2.51x | 36.9 |
| 50 (1) | 10 | 41.3 | 20.6 | 50.4 | 20.6 | 1.00 | 2.44x | 37.9 |
| 56 (1) | 5 | 41.4 | 20.8 | 50.8 | 20.8 | 1.00 | 2.44x | 36.2 |
| 60 (1) | 2 | 41.5 | 20.8 | 51.8 | 20.8 | 1.00 | 2.48x | 36.7 |
| 100 (2) | 5 | 69.3 | 72.5 | 49.8 | 52.2 | 0.72 | 0.95x | 38.6 |
| 100 (2) | 20 | 68.6 | 71.6 | 50.2 | 53.0 | 0.74 | 0.95x | 38.6 |

At one block the pruned scan now takes at most 1.00 times the time of plain Myers (same loop, so the difference is run-to-run noise); for several blocks the pruned scan is unchanged (m = 100 rows). The register-only loop also makes plain Myers faster for m <= 64, and the filters' window verification for short reads.


## Phase 3B, item 4: generalised pigeonhole (s required hits in a diagonal band of width k)

k + s disjoint seeds; a candidate needs s distinct seeds hit within a band of k diagonals; s = 1 is the classic filter. Pure filters (cost-model fallback off), mean per read; `windows` = merged candidate windows verified, `bp` = bases verified, `gen` = time to build candidates (sort and band pass). Every filter's matches are asserted equal to Myers on every read.


### E. coli 536, m = 100 (10% error), k = 10

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 72.333 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.953 | 508 | 493 | 59307 | 0.026 | 0.000 |
| filter_fa (fa, s=1) | 0.447 | 220 | 215 | 25771 | 0.011 | 0.012 |
| filter_s2 (equal, s=2) | 0.080 | 1368 | 5 | 678 | 0.068 | 0.000 |
| filter_s3 (equal, s=3) | 0.298 | 4813 | 1 | 145 | 0.273 | 0.000 |
| filter_fs2 (fa, s=2) | 0.047 | 678 | 2 | 241 | 0.033 | 0.011 |
| filter_fs3 (fa, s=3) | 0.102 | 1735 | 1 | 130 | 0.091 | 0.012 |

### E. coli 536, m = 100 (12% error), k = 12

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 73.110 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 8.436 | 4621 | 4023 | 532929 | 0.293 | 0.000 |
| filter_fa (fa, s=1) | 3.491 | 1802 | 1702 | 215592 | 0.103 | 0.016 |
| filter_s2 (equal, s=2) | 0.979 | 8147 | 231 | 29942 | 0.502 | 0.000 |
| filter_s3 (equal, s=3) | 1.433 | 21128 | 27 | 3548 | 1.379 | 0.000 |
| filter_fs2 (fa, s=2) | 0.318 | 3986 | 40 | 5187 | 0.219 | 0.012 |
| filter_fs3 (fa, s=3) | 0.523 | 8188 | 4 | 467 | 0.499 | 0.014 |

### E. coli 536, m = 100 (15% error), k = 15

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 72.623 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 41.318 | 33282 | 11915 | 2581677 | 2.580 | 0.000 |
| filter_fa (fa, s=1) | 24.992 | 14980 | 9720 | 1573613 | 1.092 | 0.026 |
| filter_s2 (equal, s=2) | 24.094 | 66872 | 7123 | 1266139 | 5.108 | 0.000 |
| filter_s3 (equal, s=3) | 12.990 | 98311 | 2474 | 368362 | 7.772 | 0.000 |
| filter_fs2 (fa, s=2) | 5.874 | 26365 | 1839 | 262385 | 1.785 | 0.017 |
| filter_fs3 (fa, s=3) | 3.885 | 44034 | 339 | 47615 | 3.077 | 0.017 |

### E. coli 536, m = 100 (17% error), k = 17

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 71.361 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 73.936 | 97109 | 7007 | 4550851 | 8.070 | 0.000 |
| filter_fa (fa, s=1) | 53.178 | 43480 | 12610 | 3303363 | 3.382 | 0.029 |
| filter_s2 (equal, s=2) | 48.011 | 120521 | 10765 | 2578861 | 9.966 | 0.000 |
| filter_s3 (equal, s=3) | 53.154 | 232017 | 4772 | 2163209 | 19.552 | 0.000 |
| filter_fs2 (fa, s=2) | 25.956 | 67827 | 7723 | 1384864 | 5.140 | 0.027 |
| filter_fs3 (fa, s=3) | 16.642 | 102493 | 3380 | 549434 | 8.163 | 0.029 |

### E. coli 536, m = 100 (20% error), k = 20

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 71.859 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 108.737 | 400572 | 521 | 4913075 | 36.346 | 0.000 |
| filter_fa (fa, s=1) | 83.893 | 150806 | 2859 | 4813202 | 12.699 | 0.034 |
| filter_s2 (equal, s=2) | 114.126 | 463463 | 222 | 4932888 | 42.069 | 0.000 |
| filter_s3 (equal, s=3) | 116.240 | 484897 | 2199 | 4837082 | 44.267 | 0.000 |
| filter_fs2 (fa, s=2) | 83.349 | 213659 | 6764 | 4387587 | 18.390 | 0.031 |
| filter_fs3 (fa, s=3) | 82.580 | 292128 | 9138 | 3761239 | 25.603 | 0.032 |

### E. coli 536, m = 100 (25% error), k = 25

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 72.922 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 244.449 | 1605638 | 1 | 4938919 | 180.121 | 0.000 |
| filter_fa (fa, s=1) | 131.846 | 633545 | 2 | 4938912 | 60.210 | 0.036 |
| filter_s2 (equal, s=2) | 284.170 | 2089283 | 2 | 4938903 | 208.657 | 0.000 |
| filter_s3 (equal, s=3) | 303.714 | 2267775 | 1 | 4938920 | 228.549 | 0.000 |
| filter_fs2 (fa, s=2) | 160.971 | 809200 | 3 | 4938880 | 79.309 | 0.038 |
| filter_fs3 (fa, s=3) | 179.661 | 1014714 | 6 | 4938806 | 102.399 | 0.040 |

### mt_orang segments vs mt_human (real homologs), m = 100, k = 15

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.245 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.173 | 133 | 45 | 10074 | 0.006 | 0.000 |
| filter_fa (fa, s=1) | 0.153 | 71 | 38 | 6528 | 0.005 | 0.019 |
| filter_s2 (equal, s=2) | 0.179 | 474 | 40 | 10018 | 0.024 | 0.000 |
| filter_s3 (equal, s=3) | 0.078 | 503 | 20 | 3345 | 0.025 | 0.000 |
| filter_fs2 (fa, s=2) | 0.041 | 123 | 10 | 1496 | 0.005 | 0.011 |
| filter_fs3 (fa, s=3) | 0.042 | 202 | 4 | 492 | 0.008 | 0.012 |

### mt_orang segments vs mt_human (real homologs), m = 100, k = 20

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.319 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.364 | 2046 | 1 | 16569 | 0.113 | 0.000 |
| filter_fa (fa, s=1) | 0.293 | 665 | 5 | 16429 | 0.044 | 0.024 |
| filter_s2 (equal, s=2) | 0.372 | 2139 | 1 | 16562 | 0.130 | 0.000 |
| filter_s3 (equal, s=3) | 0.380 | 2239 | 3 | 16492 | 0.136 | 0.000 |
| filter_fs2 (fa, s=2) | 0.298 | 931 | 13 | 15854 | 0.051 | 0.014 |
| filter_fs3 (fa, s=3) | 0.330 | 1273 | 21 | 14905 | 0.073 | 0.015 |

### mt_orang segments vs mt_human (real homologs), m = 100, k = 25

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.245 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.842 | 8985 | 1 | 16569 | 0.611 | 0.000 |
| filter_fa (fa, s=1) | 0.421 | 2748 | 1 | 16569 | 0.158 | 0.019 |
| filter_s2 (equal, s=2) | 0.945 | 9344 | 1 | 16569 | 0.820 | 0.000 |
| filter_s3 (equal, s=3) | 0.943 | 9680 | 1 | 16569 | 0.900 | 0.000 |
| filter_fs2 (fa, s=2) | 0.481 | 3438 | 1 | 16569 | 0.215 | 0.018 |
| filter_fs3 (fa, s=3) | 0.556 | 4282 | 1 | 16569 | 0.294 | 0.021 |

### mt_orang segments vs mt_human (real homologs), m = 100, k = 30

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.305 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 1.097 | 10716 | 1 | 16569 | 0.814 | 0.000 |
| filter_fa (fa, s=1) | 0.806 | 7622 | 1 | 16569 | 0.502 | 0.022 |
| filter_s2 (equal, s=2) | 1.045 | 11050 | 1 | 16569 | 0.816 | 0.000 |
| filter_s3 (equal, s=3) | 1.065 | 11385 | 1 | 16569 | 0.805 | 0.000 |
| filter_fs2 (fa, s=2) | 0.881 | 8935 | 1 | 16569 | 0.638 | 0.023 |
| filter_fs3 (fa, s=3) | 1.002 | 10351 | 1 | 16569 | 0.791 | 0.025 |

### mt_orang segments vs mt_human (real homologs), m = 200, k = 30

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.478 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.484 | 267 | 7 | 16060 | 0.012 | 0.000 |
| filter_fa (fa, s=1) | 0.409 | 106 | 21 | 12394 | 0.005 | 0.036 |
| filter_s2 (equal, s=2) | 0.352 | 276 | 20 | 11305 | 0.015 | 0.000 |
| filter_s3 (equal, s=3) | 0.149 | 284 | 13 | 4392 | 0.014 | 0.000 |
| filter_fs2 (fa, s=2) | 0.204 | 145 | 15 | 5452 | 0.007 | 0.035 |
| filter_fs3 (fa, s=3) | 0.107 | 194 | 6 | 1970 | 0.009 | 0.034 |

### mt_orang segments vs mt_human (real homologs), m = 200, k = 40

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.516 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.738 | 4094 | 1 | 16569 | 0.252 | 0.000 |
| filter_fa (fa, s=1) | 0.620 | 1150 | 1 | 16569 | 0.063 | 0.053 |
| filter_s2 (equal, s=2) | 0.780 | 4181 | 1 | 16569 | 0.274 | 0.000 |
| filter_s3 (equal, s=3) | 0.830 | 4284 | 1 | 16569 | 0.279 | 0.000 |
| filter_fs2 (fa, s=2) | 0.620 | 1369 | 1 | 16568 | 0.082 | 0.051 |
| filter_fs3 (fa, s=3) | 0.647 | 1620 | 1 | 16566 | 0.102 | 0.057 |


## Phase 3B, item 5: gated window objective (filter_fa_wg)

### Gate: fraction of reads flagged repeat-heavy (GATE_RATIO = 20, GATE_MIN = 100) and the ratio max seed count / random expectation (median, 90th percentile)

| scenario | k | reads flagged | ratio median | ratio p90 |
|---|---|---|---|---|
| E. coli uniform | 6 | 0 of 100 | 5.8 | 10.0 |
| E. coli uniform | 10 | 0 of 100 | 3.2 | 5.0 |
| E. coli uniform | 15 | 0 of 100 | 2.4 | 3.1 |
| stress uniform | 6 | 0 of 100 | 5.8 | 9.2 |
| stress uniform | 10 | 1 of 100 | 3.4 | 5.5 |
| stress uniform | 15 | 0 of 100 | 2.5 | 4.6 |
| stress clustered | 6 | 1 of 100 | 5.8 | 9.2 |
| stress clustered | 10 | 1 of 100 | 3.4 | 6.3 |
| stress clustered | 15 | 0 of 100 | 2.4 | 4.2 |
| stress low-complexity edge | 6 | 97 of 100 | 4843.3 | 5655.8 |
| stress low-complexity edge | 10 | 94 of 100 | 189.0 | 362.0 |
| stress low-complexity edge | 15 | 0 of 100 | 5.3 | 6.4 |
| stress low-complexity inside | 6 | 100 of 100 | 3931.7 | 5655.8 |
| stress low-complexity inside | 10 | 96 of 100 | 188.2 | 362.0 |
| stress low-complexity inside | 15 | 0 of 100 | 5.4 | 6.4 |

### Before and after (pure filters, mean ms per read, best of 3; selection ms in brackets)

| scenario | k | filter (equal split) | filter_fa (hits) | filter_fa_w, sort estimate (phase 2) | filter_fa_w, binned estimate | filter_fa_wg (gated, binned) | filter_fa_wg (gated, sort estimate) | filter_auto |
|---|---|---|---|---|---|---|---|---|
| stress low-complexity edge | 6 | 0.478 (0.000) | 0.042 (0.005) | 11.071 (11.034) | 0.679 (0.644) | 0.655 (0.622) | 11.974 (11.929) | 0.040 (0.011) |
| stress low-complexity edge | 10 | 1.487 (0.000) | 2.551 (0.014) | 9.538 (8.640) | 1.452 (0.529) | 1.462 (0.519) | 9.602 (8.663) | 0.444 (0.025) |
| stress low-complexity edge | 15 | 32.910 (0.000) | 25.433 (0.025) | 32.193 (4.386) | 29.158 (0.362) | 24.973 (0.022) | 24.769 (0.021) | 6.552 (0.037) |
| stress low-complexity inside | 6 | 0.781 (0.000) | 0.572 (0.005) | 23.633 (23.058) | 1.824 (1.273) | 1.848 (1.301) | 23.658 (23.070) | 0.595 (0.012) |
| stress low-complexity inside | 10 | 1.784 (0.000) | 1.565 (0.010) | 16.845 (15.525) | 2.233 (0.922) | 2.306 (0.931) | 17.071 (15.706) | 1.166 (0.024) |
| stress low-complexity inside | 15 | 25.520 (0.000) | 19.708 (0.019) | 32.911 (9.011) | 24.778 (0.592) | 19.500 (0.017) | 19.524 (0.018) | 9.654 (0.036) |
| E. coli uniform | 6 | 0.032 (0.000) | 0.013 (0.006) | 0.129 (0.122) | 0.051 (0.044) | 0.013 (0.006) | 0.014 (0.006) | 0.013 (0.008) |
| E. coli uniform | 10 | 0.937 (0.000) | 0.424 (0.012) | 0.544 (0.128) | 0.481 (0.066) | 0.428 (0.013) | 0.423 (0.012) | 0.065 (0.024) |
| E. coli uniform | 15 | 39.430 (0.000) | 24.856 (0.026) | 24.779 (0.190) | 25.450 (0.120) | 25.065 (0.026) | 24.780 (0.025) | 3.861 (0.038) |
| stress clustered | 6 | 0.046 (0.000) | 0.023 (0.006) | 0.345 (0.328) | 0.080 (0.062) | 0.038 (0.021) | 0.224 (0.207) | 0.026 (0.012) |
| stress clustered | 10 | 1.095 (0.000) | 0.504 (0.013) | 1.199 (0.690) | 0.588 (0.091) | 0.529 (0.038) | 1.040 (0.546) | 0.097 (0.024) |
| stress clustered | 15 | 38.328 (0.000) | 25.669 (0.030) | 26.831 (0.833) | 26.673 (0.161) | 26.992 (0.028) | 26.290 (0.028) | 4.202 (0.039) |


## Phase 3B: s-hit filters inside filter_auto

Mean ms per read (fallback active), best of 3 runs, 100 reads. Choices: queries that filter_auto sent to the 2-hit / 3-hit filter (of all queries).

| text | setting | Myers | filter | filter_fa | filter_fs2 | filter_fs3 | filter_auto (no s) | filter_auto | auto / best fixed | choices s2 / s3 |
|---|---|---|---|---|---|---|---|---|---|---|
| E. coli 536 | m=100, k=6 | 74.765 | 0.033 | 0.013 | 0.010 | 0.010 | 0.013 | 0.014 | 1.43 | 22 / 0 |
| E. coli 536 | m=100, k=8 | 77.427 | 0.059 | 0.034 | 0.013 | 0.023 | 0.034 | 0.026 | 1.91 | 98 / 0 |
| E. coli 536 | m=100, k=10 | 73.235 | 0.928 | 0.432 | 0.046 | 0.104 | 0.445 | 0.065 | 1.40 | 100 / 0 |
| E. coli 536 | m=100, k=12 | 74.776 | 8.556 | 3.390 | 0.333 | 0.616 | 3.405 | 0.375 | 1.13 | 100 / 0 |
| E. coli 536 | m=100, k=14 | 73.575 | 31.201 | 14.142 | 2.158 | 1.856 | 13.928 | 2.041 | 1.10 | 57 / 43 |
| E. coli 536 | m=100, k=15 | 73.273 | 38.527 | 25.231 | 5.628 | 3.738 | 24.544 | 3.956 | 1.06 | 2 / 98 |
| E. coli 536 | m=100, k=17 | 73.675 | 74.167 | 55.278 | 26.708 | 16.429 | 53.983 | 17.031 | 1.04 | 0 / 100 |
| E. coli 536 | m=100, k=20 | 75.518 | 74.058 | 76.667 | 90.208 | 92.481 | 75.226 | 73.844 | 1.00 | 0 / 7 |
| E. coli 536 | m=200, k=20 | 145.671 | 4.305 | 1.598 | 0.065 | 0.077 | 1.616 | 0.114 | 1.74 | 100 / 0 |
| E. coli 536 | m=200, k=30 | 144.487 | 129.787 | 97.209 | 31.314 | 10.210 | 97.877 | 10.055 | 0.98 | 0 / 100 |
| E. coli 536 | m=400, k=60 | 261.226 | 266.281 | 250.991 | 194.847 | 101.960 | 253.869 | 101.470 | 1.00 | 0 / 97 |
| synthetic | m=100, k=10 | 67.445 | 0.634 | 0.357 | 0.042 | 0.096 | 0.372 | 0.061 | 1.48 | 100 / 0 |
| synthetic | m=100, k=12 | 67.608 | 6.279 | 3.169 | 0.276 | 0.496 | 3.170 | 0.310 | 1.12 | 100 / 0 |
| synthetic | m=100, k=15 | 68.191 | 33.563 | 25.359 | 6.328 | 4.350 | 25.144 | 4.481 | 1.03 | 1 / 99 |
| mt_orang vs mt_human (real) | m=100, k=15 | 0.248 | 0.159 | 0.117 | 0.041 | 0.029 | 0.115 | 0.052 | 1.78 | 0 / 100 |
| mt_orang vs mt_human (real) | m=100, k=20 | 0.253 | 0.243 | 0.258 | 0.312 | 0.331 | 0.266 | 0.277 | 1.14 | 0 / 0 |
| mt_orang vs mt_human (real) | m=200, k=30 | 0.482 | 0.491 | 0.411 | 0.205 | 0.104 | 0.413 | 0.173 | 1.66 | 0 / 100 |

Worst auto/best-fixed ratio: 1.91 (E. coli 536, m=100, k=8); 12 of 17 settings exceed 1.05.


## Final pass: interleaved per-read timing, speedup over Myers

Every method runs on every read inside one process, in a rotating order; 3 repeats of 60 reads per setting; index build and process startup are excluded. Cell = median speedup over Myers (time(Myers)/time(method) on the same read), with the 10th to 90th percentile in brackets (>1 faster than the scan). All methods return exactly the same end positions in every run (checked per read). Fallback active. Parameters frozen. The blind table is the B. fragilis one.


### E. coli 536 (tuned on)

| m, k | lv | filter_fa | filter_fa_wg | filter_s2 | filter_s3 | filter_fs2 | filter_fs3 | filter_auto |
|---|---|---|---|---|---|---|---|---|
| m=100, k=1 | 4.77 (4.19 to 5.08) | 3255.87 (2731.87 to 3913.64) | 14528.17 (3393.35 to 17415.52) | 18739.68 (14053.15 to 23584.50) | 28183.59 (17484.64 to 33262.26) | 14982.14 (4831.52 to 16945.92) | 15299.17 (5593.87 to 17269.58) | 14912.26 (6632.59 to 18412.06) |
| m=100, k=2 | 3.29 (2.91 to 3.48) | 3119.89 (2589.67 to 3558.50) | 13157.36 (3217.36 to 15120.04) | 17476.10 (13050.25 to 21811.27) | 23533.13 (15144.53 to 30060.69) | 13612.66 (4833.26 to 16064.11) | 14134.54 (5296.19 to 16236.28) | 13006.41 (6180.63 to 15811.20) |
| m=100, k=4 | 1.97 (1.68 to 2.11) | 2884.35 (2324.68 to 3334.61) | 9957.87 (3054.57 to 11613.91) | 16785.53 (13211.45 to 20915.73) | 21289.77 (14466.09 to 27001.20) | 12234.35 (4642.25 to 14124.55) | 12796.68 (5027.48 to 14592.29) | 10566.01 (5520.52 to 13097.95) |
| m=100, k=6 | 1.44 (1.31 to 1.52) | 2442.39 (1884.90 to 2913.52) | 7408.80 (2640.53 to 10036.40) | 14722.60 (11283.06 to 18989.73) | 17680.79 (9856.09 to 24040.10) | 10873.16 (4275.90 to 12588.83) | 10315.45 (4851.62 to 12034.67) | 7900.01 (3517.72 to 10995.95) |
| m=100, k=8 | 1.12 (1.03 to 1.18) | 1357.61 (946.17 to 1825.37) | 2601.81 (1282.71 to 4199.39) | 5524.92 (3229.52 to 9869.63) | 3102.63 (1025.97 to 4664.66) | 3960.20 (2906.60 to 5132.30) | 2931.10 (2153.23 to 3581.63) | 3293.55 (1781.88 to 4001.68) |
| m=100, k=10 | 0.92 (0.85 to 0.97) | 156.53 (95.65 to 230.47) | 178.29 (105.35 to 284.91) | 918.80 (616.46 to 1277.74) | 237.62 (188.03 to 305.33) | 1414.22 (984.85 to 2001.90) | 653.38 (433.78 to 940.17) | 1100.00 (768.10 to 1461.80) |
| m=100, k=12 | 0.78 (0.72 to 0.82) | 19.83 (14.24 to 32.48) | 20.96 (14.91 to 34.83) | 141.68 (26.26 to 193.68) | 46.77 (37.39 to 58.33) | 214.90 (129.70 to 393.63) | 134.44 (90.47 to 204.84) | 205.94 (127.37 to 355.98) |
| m=100, k=15 | 0.64 (0.56 to 0.68) | 2.87 (2.13 to 3.92) | 2.89 (2.08 to 4.03) | 2.12 (1.73 to 16.10) | 5.61 (4.46 to 6.79) | 13.28 (7.12 to 23.89) | 19.49 (11.26 to 30.04) | 18.99 (11.16 to 29.33) |
| m=100, k=17 | 0.58 (0.52 to 0.63) | 1.37 (1.15 to 1.73) | 1.35 (1.14 to 1.75) | 1.60 (1.32 to 2.03) | 3.42 (0.64 to 4.65) | 2.96 (1.84 to 4.98) | 5.03 (2.54 to 8.14) | 5.08 (2.54 to 8.69) |
| m=100, k=20 | 0.49 (0.42 to 0.51) | 1.00 (0.93 to 1.06) | 1.00 (0.90 to 1.05) | 0.64 (0.57 to 0.68) | 0.62 (0.56 to 0.68) | 0.82 (0.74 to 1.05) | 0.93 (0.69 to 1.22) | 0.99 (0.93 to 1.15) |
| m=100, k=25 | 0.40 (0.36 to 0.42) | 1.00 (0.91 to 1.07) | 1.00 (0.95 to 1.07) | 1.00 (0.86 to 1.07) | 1.00 (0.89 to 1.07) | 0.48 (0.38 to 0.57) | 0.43 (0.35 to 0.59) | 1.00 (0.95 to 1.04) |
| m=200, k=30 | 0.64 (0.59 to 0.73) | 1.51 (1.28 to 1.89) | 1.52 (1.25 to 1.84) | 2.31 (1.98 to 2.94) | 8.57 (1.07 to 12.15) | 4.67 (2.99 to 8.33) | 16.24 (8.86 to 33.47) | 16.12 (8.93 to 32.72) |
| m=200, k=40 | 0.50 (0.43 to 0.53) | 1.00 (0.94 to 1.07) | 1.00 (0.94 to 1.07) | 0.63 (0.57 to 0.69) | 0.62 (0.56 to 0.68) | 0.84 (0.75 to 0.91) | 0.81 (0.72 to 0.88) | 0.99 (0.95 to 1.03) |

Median absolute time of the Myers scan (ms per read, for scale only): m=100,k=1: 73.3, m=100,k=2: 73.1, m=100,k=4: 73.8, m=100,k=6: 73.0, m=100,k=8: 72.5, m=100,k=10: 72.6, m=100,k=12: 72.4, m=100,k=15: 73.1, m=100,k=17: 74.2, m=100,k=20: 72.6, m=100,k=25: 72.8, m=200,k=30: 143.2, m=200,k=40: 142.0.


### B. fragilis (BLIND)

| m, k | lv | filter_fa | filter_fa_wg | filter_s2 | filter_s3 | filter_fs2 | filter_fs3 | filter_auto |
|---|---|---|---|---|---|---|---|---|
| m=100, k=1 | 3.87 (3.63 to 4.06) | 3425.61 (2935.45 to 4058.76) | 15632.46 (3582.99 to 18090.93) | 18940.79 (14951.73 to 25672.48) | 29792.62 (17062.28 to 35255.43) | 15676.60 (5418.81 to 17720.12) | 15964.03 (5901.28 to 17771.80) | 15671.21 (8833.19 to 18752.85) |
| m=100, k=2 | 3.22 (2.60 to 3.45) | 3255.75 (2696.44 to 3674.40) | 13004.10 (3160.73 to 15592.98) | 19405.79 (14935.31 to 23564.15) | 25628.76 (16758.20 to 31823.11) | 14300.31 (5151.08 to 16721.79) | 14616.97 (5660.86 to 16523.94) | 13739.15 (6275.14 to 15591.22) |
| m=100, k=4 | 1.97 (1.79 to 2.09) | 2899.00 (2311.77 to 3260.95) | 10824.71 (2982.79 to 12603.41) | 16614.29 (12440.47 to 20445.19) | 21994.04 (15082.92 to 27504.10) | 12838.53 (4763.17 to 14247.15) | 13281.28 (5000.58 to 14874.62) | 11647.98 (6324.37 to 13190.82) |
| m=100, k=6 | 1.41 (1.25 to 1.52) | 2471.41 (1907.23 to 2886.14) | 7385.06 (2626.91 to 9054.73) | 15839.46 (11924.67 to 19780.47) | 19409.78 (10133.08 to 24268.84) | 10320.82 (4568.66 to 12989.20) | 10216.91 (4778.08 to 12671.51) | 7637.50 (3882.82 to 10310.60) |
| m=100, k=8 | 1.15 (0.98 to 1.21) | 1277.82 (793.04 to 1767.88) | 2357.87 (1080.42 to 4180.39) | 6559.13 (3203.73 to 10742.40) | 3593.68 (1082.40 to 5404.12) | 3899.73 (2680.48 to 5331.92) | 2988.32 (1893.22 to 3894.99) | 3308.86 (1906.36 to 4035.70) |
| m=100, k=10 | 0.91 (0.79 to 0.98) | 146.16 (90.34 to 258.03) | 167.01 (100.89 to 307.66) | 914.59 (517.18 to 1186.30) | 241.57 (184.13 to 387.71) | 1386.02 (842.01 to 2027.31) | 645.66 (389.99 to 1019.85) | 1030.38 (638.76 to 1582.23) |
| m=100, k=12 | 0.77 (0.68 to 0.83) | 20.78 (13.86 to 30.50) | 22.26 (14.51 to 31.97) | 147.26 (27.74 to 203.79) | 49.72 (36.18 to 65.00) | 224.62 (139.39 to 326.53) | 134.35 (90.38 to 196.87) | 207.80 (127.31 to 308.27) |
| m=100, k=15 | 0.63 (0.58 to 0.67) | 2.74 (2.17 to 4.12) | 2.72 (2.15 to 4.18) | 2.13 (1.61 to 14.05) | 5.32 (3.72 to 6.66) | 11.12 (7.64 to 21.62) | 17.30 (11.93 to 32.01) | 17.17 (12.16 to 29.14) |
| m=100, k=17 | 0.57 (0.51 to 0.61) | 1.31 (1.07 to 1.60) | 1.33 (1.08 to 1.61) | 1.59 (1.24 to 2.07) | 2.67 (0.63 to 4.64) | 2.65 (1.72 to 4.32) | 4.26 (2.49 to 7.72) | 4.20 (2.44 to 7.71) |
| m=100, k=20 | 0.47 (0.41 to 0.51) | 1.00 (0.90 to 1.12) | 1.00 (0.90 to 1.11) | 0.63 (0.57 to 0.72) | 0.62 (0.55 to 0.70) | 0.79 (0.69 to 1.01) | 0.75 (0.63 to 1.16) | 0.99 (0.94 to 1.16) |
| m=100, k=25 | 0.39 (0.34 to 0.43) | 1.00 (0.89 to 1.13) | 1.00 (0.92 to 1.14) | 1.00 (0.87 to 1.12) | 1.00 (0.87 to 1.13) | 0.46 (0.38 to 0.54) | 0.41 (0.34 to 0.50) | 1.00 (0.94 to 1.06) |
| m=200, k=30 | 0.65 (0.54 to 0.71) | 1.45 (1.23 to 1.78) | 1.44 (1.22 to 1.74) | 2.26 (1.79 to 2.77) | 8.49 (1.00 to 13.63) | 4.62 (2.69 to 6.23) | 15.18 (7.31 to 23.77) | 15.07 (7.47 to 22.16) |
| m=200, k=40 | 0.49 (0.40 to 0.54) | 1.00 (0.90 to 1.14) | 1.00 (0.88 to 1.15) | 0.63 (0.55 to 0.81) | 0.62 (0.54 to 0.71) | 0.84 (0.72 to 0.97) | 0.81 (0.72 to 0.93) | 0.99 (0.91 to 1.07) |

Median absolute time of the Myers scan (ms per read, for scale only): m=100,k=1: 76.0, m=100,k=2: 76.7, m=100,k=4: 75.3, m=100,k=6: 76.0, m=100,k=8: 76.1, m=100,k=10: 77.3, m=100,k=12: 75.9, m=100,k=15: 76.3, m=100,k=17: 75.8, m=100,k=20: 76.6, m=100,k=25: 76.0, m=200,k=30: 149.8, m=200,k=40: 150.0.



## Final pass: gate on B. fragilis (BLIND, thresholds frozen: ratio 20, minimum count 100)

| k | reads flagged | ratio median | ratio p90 | ratio max |
|---|---|---|---|---|
| 6 | 0 of 200 | 5.0 | 9.2 | 32.5 |
| 10 | 0 of 200 | 3.3 | 5.4 | 13.0 |
| 15 | 0 of 200 | 2.7 | 3.9 | 4.5 |

