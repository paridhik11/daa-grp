# Benchmark summary

Main text: **E. coli 536** (see `data/SOURCES.txt`; no K-12 genome was reachable). Search time is per pattern, single thread; index build is reported separately.

Edlib and WFA2 return the BEST match per read (minimum distance, ends at that minimum), not all matches within k, and include Python call overhead; compare them on time only.

Cost-model constants (ns), calibrated on synthetic genome_synth.txt and evaluated on E. coli 536: c_hit=109.7, c_win=219, c_bp=7.28, a_scan=7.70, c_sel=57.7


Headline experiment (E. coli 536, m = 100, 100 reads per setting): values are the **median [10th, 90th percentile] over reads** of per-read query time; other experiments use means over 20 reads (Sellers 1, WFA2 5).


## Frequency-aware vs equal split (pure filter, no fallback), E. coli 536

| k | hits equal | hits fa | hit ratio (medians) | bp verified equal | bp verified fa | time equal (ms) | time fa (ms) | fa speedup vs equal | fa speedup vs Myers |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 5 [2, 10] | 1 [1, 2] | 5.0x | 410 [104, 927] | 102 [101, 103] | 0.009 [0.003, 0.022] | 0.028 [0.025, 0.034] | 0.33x | 2563x |
| 2 | 9 [5, 16] | 1 [1, 2] | 9.0x | 736 [404, 1477] | 104 [103, 208] | 0.013 [0.007, 0.028] | 0.030 [0.028, 0.040] | 0.43x | 2409x |
| 3 | 12 [6, 18] | 2 [1, 3] | 6.0x | 1041 [526, 1739] | 107 [105, 214] | 0.024 [0.012, 0.043] | 0.036 [0.032, 0.044] | 0.67x | 1981x |
| 4 | 14 [9, 25] | 3 [1, 4] | 4.5x | 1242 [641, 2424] | 110 [107, 322] | 0.021 [0.012, 0.040] | 0.037 [0.034, 0.044] | 0.57x | 1937x |
| 6 | 18 [11, 31] | 4 [3, 7] | 4.5x | 1632 [1008, 3145] | 237 [114, 666] | 0.030 [0.018, 0.055] | 0.045 [0.041, 0.056] | 0.68x | 1550x |
| 8 | 24 [14, 68] | 15 [8, 25] | 1.6x | 2437 [1381, 6990] | 1400 [708, 2439] | 0.053 [0.031, 0.163] | 0.067 [0.054, 0.087] | 0.79x | 1053x |
| 10 | 326 [216, 1184] | 212 [115, 321] | 1.5x | 38940 [26171, 134831] | 24740 [13387, 37104] | 0.675 [0.469, 2.679] | 0.466 [0.296, 0.713] | 1.45x | 153x |
| 11 | 1260 [978, 1592] | 667 [347, 1068] | 1.9x | 150050 [116381, 187210] | 80439 [43311, 125188] | 2.583 [1.978, 3.346] | 1.759 [0.986, 2.775] | 1.47x | 41x |
| 12 | 4828 [1813, 6220] | 1656 [1071, 2639] | 2.9x | 558434 [227061, 711286] | 204320 [133726, 310830] | 8.614 [3.440, 11.021] | 3.258 [2.149, 4.925] | 2.64x | 22x |
| 13 | 5494 [4530, 19394] | 3774 [2303, 6031] | 1.5x | 654841 [544943, 1873433] | 456108 [288378, 689800] | 10.775 [8.788, 29.522] | 7.067 [4.550, 11.241] | 1.52x | 10x |
| 14 | 20796 [13438, 26218] | 7370 [5305, 11344] | 2.8x | 2041166 [1413272, 2405927] | 867459 [639550, 1244160] | 31.989 [22.387, 39.469] | 13.384 [9.799, 18.451] | 2.39x | 5x |
| 15 | 23647 [19932, 28911] | 14000 [9854, 20383] | 1.7x | 2277172 [2022391, 2639100] | 1523332 [1110920, 1984815] | 35.187 [31.098, 43.073] | 23.837 [17.912, 31.487] | 1.48x | 3x |
| 17 | 97332 [83921, 111570] | 41428 [28912, 59795] | 2.3x | 4570596 [4396164, 4663902] | 3291704 [2696102, 3903638] | 74.113 [69.708, 83.170] | 53.071 [42.207, 67.460] | 1.40x | 1x |
| 20 | 435776 [128342, 492056] | 144428 [108145, 209961] | 3.0x | 4938868 [4816463, 4938916] | 4849136 [4684707, 4914649] | 113.957 [81.705, 123.068] | 83.613 [77.223, 94.813] | 1.36x | 1x |
| 30 | 2504724 [2330961, 2638079] | 1747860 [1409690, 2349176] | 1.4x | 4938920 [4938920, 4938920] | 4938920 [4938920, 4938920] | 360.144 [336.338, 401.907] | 264.109 [219.684, 341.831] | 1.36x | 0x |

## Cost model (predicted vs measured), E. coli 536

Fit: NNLS (scipy.optimize.nnls) on synthetic genome_synth.txt; evaluated here on E. coli 536 (a different text, but NOT blind: the calibration range was widened after a first fit failed on it; see the blind-evaluation section). c_win = 219 ns (unconstrained least squares gave 219, so the sign constraint was not active). Hit counts come from the index; windows and merged bp are estimated assuming uniform hit placement.

- Crossover filter vs Myers: measured k/m = 0.168 (k = 16.8), predicted k/m = 0.169 (k = 16.9)
- Crossover filter_fa vs Myers: measured k/m = 0.187 (k = 18.7), predicted k/m = 0.191 (k = 19.1)
- filter time prediction: median factor error 1.09x, worst 1.27x over 14 values of k (measured median > 0.01 ms; prediction ratio uses medians)
- filter_fa time prediction: median factor error 1.12x, worst 2.36x over 15 values of k (measured median > 0.01 ms; prediction ratio uses medians)
- Myers scan prediction: median factor error 1.07x, worst 1.11x

| k | filter measured (ms) | predicted | filter_fa measured | predicted | Myers measured | predicted |
|---|---|---|---|---|---|---|
| 1 | 0.009 [0.003, 0.022] | 0.009 | 0.028 [0.025, 0.034] | 0.065 | 70.8 [67.6, 83.5] | 76.0 |
| 2 | 0.013 [0.007, 0.028] | 0.016 | 0.030 [0.028, 0.040] | 0.067 | 73.0 [67.9, 87.4] | 76.0 |
| 3 | 0.024 [0.012, 0.043] | 0.022 | 0.036 [0.032, 0.044] | 0.067 | 70.8 [67.8, 86.6] | 76.0 |
| 4 | 0.021 [0.012, 0.040] | 0.026 | 0.037 [0.034, 0.044] | 0.069 | 71.4 [67.8, 82.0] | 76.0 |
| 6 | 0.030 [0.018, 0.055] | 0.036 | 0.045 [0.041, 0.056] | 0.072 | 69.7 [65.7, 80.2] | 76.0 |
| 8 | 0.053 [0.031, 0.163] | 0.048 | 0.067 [0.054, 0.087] | 0.094 | 71.0 [67.8, 77.4] | 76.0 |
| 10 | 0.675 [0.469, 2.679] | 0.680 | 0.466 [0.296, 0.713] | 0.500 | 71.3 [67.7, 88.6] | 76.0 |
| 11 | 2.583 [1.978, 3.346] | 2.597 | 1.759 [0.986, 2.775] | 1.456 | 72.4 [69.0, 87.5] | 76.0 |
| 12 | 8.614 [3.440, 11.021] | 9.627 | 3.258 [2.149, 4.925] | 3.578 | 70.2 [66.9, 77.6] | 76.0 |
| 13 | 10.775 [8.788, 29.522] | 11.210 | 7.067 [4.550, 11.241] | 7.893 | 70.9 [68.2, 89.3] | 76.0 |
| 14 | 31.989 [22.387, 39.469] | 34.810 | 13.384 [9.799, 18.451] | 14.794 | 69.3 [66.2, 75.5] | 76.0 |
| 15 | 35.187 [31.098, 43.073] | 38.912 | 23.837 [17.912, 31.487] | 26.136 | 69.6 [66.5, 74.1] | 76.0 |
| 17 | 74.113 [69.708, 83.170] | 79.096 | 53.071 [42.207, 67.460] | 56.331 | 68.7 [65.6, 77.5] | 76.0 |
| 20 | 113.957 [81.705, 123.068] | 119.680 | 83.613 [77.223, 94.813] | 87.127 | 68.9 [66.5, 71.9] | 76.0 |
| 30 | 360.144 [336.338, 401.907] | 346.593 | 264.109 [219.684, 341.831] | 263.647 | 70.3 [66.8, 77.7] | 76.0 |

## Index build vs query time and break-even, E. coli 536

Q* = index_ms / (Myers median ms per query - filter median ms per query); the number of queries after which building the index plus filtering beats repeated plain Myers scans (never, if the filter is not faster). Query times are medians over 100 reads (with the cost-model fallback active); the index build time is one measurement per setting.

| k | index build (ms) | Myers (ms/query) | filter (ms/query) | Q* filter | filter_fa (ms/query) | Q* filter_fa |
|---|---|---|---|---|---|---|
| 1 | 110 | 70.8 | 0.007 | 1.55 | 0.029 | 1.76 |
| 2 | 142 | 73.0 | 0.017 | 1.95 | 0.034 | 1.49 |
| 3 | 148 | 70.8 | 0.023 | 2.09 | 0.034 | 1.86 |
| 4 | 106 | 71.4 | 0.021 | 1.48 | 0.036 | 1.47 |
| 6 | 144 | 69.7 | 0.029 | 2.06 | 0.043 | 1.46 |
| 8 | 127 | 71.0 | 0.040 | 1.78 | 0.084 | 1.95 |
| 10 | 116 | 71.3 | 0.631 | 1.65 | 0.455 | 1.81 |
| 11 | 118 | 72.4 | 2.581 | 1.70 | 1.494 | 1.66 |
| 12 | 114 | 70.2 | 9.117 | 1.86 | 3.238 | 1.65 |
| 13 | 117 | 70.9 | 10.709 | 1.95 | 7.256 | 1.77 |
| 14 | 120 | 69.3 | 32.443 | 3.26 | 13.982 | 1.92 |
| 15 | 111 | 69.6 | 35.638 | 3.26 | 23.340 | 2.28 |
| 17 | 107 | 68.7 | 70.546 | inf | 53.992 | 8.15 |
| 20 | 107 | 68.9 | 71.311 | inf | 72.574 | inf |
| 30 | 105 | 70.3 | 70.945 | inf | 71.556 | inf |

## Cost model: which figures are blind and which are fitted

Protocol (fixed before the blind text was read): constants are fitted by NNLS on the synthetic genome X (synthetic genome_synth.txt) with k up to 25, then frozen. Prediction error = median over k of the factor max(pred/measured, measured/pred) between the median predicted and the median measured per-read time of the pure filter (fallback off), for k where the measured median exceeds 0.01 ms. The fit was not changed after seeing any evaluation text.

| text | status | reads per k | filter error (median / worst) | filter_fa error (median / worst) | Myers error (median) | crossover filter measured / predicted (k/m) | crossover filter_fa measured / predicted (k/m) |
|---|---|---|---|---|---|---|---|
| Synthetic genome X | fitted text, new reads | 30 | 1.13x / 1.36x | 1.10x / 1.90x | 1.08x | 0.169 / 0.169 | 0.183 / 0.185 |
| E. coli 536 | not fitted, range chosen after seeing it | 100 | 1.09x / 1.27x | 1.12x / 2.36x | 1.07x | 0.168 / 0.169 | 0.187 / 0.191 |
| B. fragilis | blind | 100 | 1.11x / 1.25x | 1.15x / 2.36x | 1.11x | 0.168 / 0.169 | 0.186 / 0.190 |

Fitted: the constants themselves and the synthetic-text row. Not blind: the E. coli row (the calibration range was widened to k <= 25 after a first fit restricted to k <= 15 failed on this text). Blind: the B. fragilis row only.


## Cells updated per second, E. coli 536 (m = 100)

Sellers and Myers: m*n DP cells represented per pattern (Myers updates 64 per word operation); banded: cells actually computed.

| k | Sellers cells/s | banded cells computed | banded cells/s | Myers cells represented | Myers cells/s | Myers/banded cells/s | banded ms | Myers ms |
|---|---|---|---|---|---|---|---|---|
| 1 | 5.38e+08 | 3.09e+08 | 2.03e+08 | 9.78e+09 | 6.16e+09 | 30x | 76.4 | 79.3 |
| 2 | 5.24e+08 | 4.91e+08 | 2.62e+08 | 9.88e+09 | 6.67e+09 | 25x | 93.7 | 74.1 |
| 3 | 5.28e+08 | 6.79e+08 | 2.82e+08 | 9.88e+09 | 6.63e+09 | 24x | 120.3 | 74.4 |
| 4 | 5.24e+08 | 8.64e+08 | 3.25e+08 | 9.78e+09 | 6.48e+09 | 20x | 133.0 | 75.5 |
| 6 | 5.23e+08 | 1.23e+09 | 3.49e+08 | 9.88e+09 | 6.85e+09 | 20x | 176.8 | 72.1 |
| 8 | 5.73e+08 | 1.66e+09 | 4.00e+08 | 9.68e+09 | 6.04e+09 | 15x | 207.5 | 80.1 |
| 10 | 5.67e+08 | 2.00e+09 | 4.14e+08 | 9.68e+09 | 6.43e+09 | 16x | 241.8 | 75.2 |
| 11 | 5.56e+08 | 2.21e+09 | 4.41e+08 | 9.58e+09 | 6.78e+09 | 15x | 250.4 | 70.7 |
| 12 | 5.70e+08 | 2.40e+09 | 4.41e+08 | 9.48e+09 | 6.61e+09 | 15x | 272.5 | 71.8 |
| 13 | 5.10e+08 | 2.60e+09 | 4.40e+08 | 9.48e+09 | 6.43e+09 | 15x | 295.5 | 73.8 |
| 14 | 5.42e+08 | 2.78e+09 | 4.43e+08 | 9.48e+09 | 6.14e+09 | 14x | 314.2 | 77.2 |
| 15 | 5.40e+08 | 2.98e+09 | 4.41e+08 | 9.58e+09 | 6.08e+09 | 14x | 337.3 | 78.8 |
| 17 | 5.50e+08 | 3.37e+09 | 4.63e+08 | 9.38e+09 | 6.49e+09 | 14x | 364.8 | 72.3 |
| 20 | 5.57e+08 | 3.97e+09 | 4.80e+08 | 9.48e+09 | 6.70e+09 | 14x | 413.4 | 70.7 |
| 30 | 5.57e+08 | 5.83e+09 | 4.91e+08 | 9.88e+09 | 6.73e+09 | 14x | 593.4 | 73.3 |

## Stress scenarios, E. coli 536 + ~1% injected microsatellite/homopolymer arrays (m = 100, 20 reads each)

Recall by read = fraction of reads (that have a true match within k) for which the heuristic reports at least one match; recall by end = fraction of true end positions reported. Exact methods are asserted equal to Myers in every row.

| scenario | k | one-hit read | one-hit end | two-hit read | two-hit end | filter hits | filter_fa hits | filter ms | filter_fa ms | Myers ms |
|---|---|---|---|---|---|---|---|---|---|---|
| uniform | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17 | 4 | 0.030 | 0.051 | 72.9 |
| uniform | 10 | 1.00 | 1.00 | 0.95 | 0.97 | 463 | 186 | 0.870 | 0.410 | 70.1 |
| uniform | 15 | 1.00 | 1.00 | 0.45 | 1.00 | 25050 | 15693 | 32.962 | 23.103 | 72.0 |
| clustered | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 22 | 8 | 0.032 | 0.048 | 67.3 |
| clustered | 10 | 1.00 | 1.00 | 1.00 | 1.00 | 3256 | 2446 | 1.211 | 0.911 | 67.9 |
| clustered | 15 | 1.00 | 1.00 | 1.00 | 1.00 | 28859 | 18291 | 35.408 | 25.864 | 71.8 |
| indel_heavy | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17 | 5 | 0.033 | 0.053 | 69.6 |
| indel_heavy | 10 | 1.00 | 1.00 | 0.90 | 0.98 | 679 | 248 | 1.307 | 0.517 | 70.5 |
| indel_heavy | 15 | 0.95 | 0.97 | 0.60 | 0.79 | 33642 | 16342 | 44.050 | 30.287 | 76.9 |
| lowcomplex_edge | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 6878 | 12 | 0.627 | 0.054 | 63.6 |
| lowcomplex_edge | 10 | 1.00 | 1.00 | 0.95 | 0.97 | 12693 | 2839 | 1.631 | 3.249 | 70.7 |
| lowcomplex_edge | 15 | 1.00 | 1.00 | 0.65 | 0.59 | 38737 | 26362 | 35.723 | 27.991 | 70.0 |
| lowcomplex_inside | 6 | 1.00 | 1.00 | 1.00 | 1.00 | 17734 | 7111 | 1.069 | 0.767 | 62.3 |
| lowcomplex_inside | 10 | 1.00 | 1.00 | 0.95 | 0.99 | 18046 | 11129 | 2.593 | 1.656 | 63.8 |
| lowcomplex_inside | 15 | 1.00 | 0.91 | 0.75 | 0.68 | 51356 | 33393 | 24.282 | 21.027 | 72.9 |

## Results by dataset (m = 100; ms per pattern; index build once)

| dataset | length (bp) | k | Sellers | Myers | Edlib (best/read) | WFA2 (best/read) | filter | filter_fa | one-hit recall (read) | two-hit recall (read) | index build (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| E. coli 536 | 4938920 | 5 | 836.05 | 75.63 | 37.43 | 183.71 | 0.0264 | 0.0465 | 1.00 | 1.00 | 101.0 |
| E. coli 536 | 4938920 | 10 | 904.53 | 74.34 | 36.54 | 322.52 | 0.7611 | 0.4012 | 1.00 | 0.95 | 108.0 |
| SARS-CoV-2 (MT192765.1) | 29829 | 5 | 5.22 | 0.51 | 0.27 | 0.82 | 0.0033 | 0.0335 | 1.00 | 1.00 | 21.3 |
| SARS-CoV-2 (MT192765.1) | 29829 | 10 | 6.17 | 0.41 | 0.24 | 1.26 | 0.0164 | 0.0455 | 1.00 | 0.95 | 24.6 |
| phage lambda | 48502 | 5 | 8.26 | 0.68 | 0.36 | 1.11 | 0.0036 | 0.0339 | 1.00 | 1.00 | 21.2 |
| phage lambda | 48502 | 10 | 8.68 | 0.69 | 0.35 | 4.19 | 0.0192 | 0.0600 | 1.00 | 0.95 | 28.0 |
| human mtDNA | 16569 | 5 | 3.22 | 0.24 | 0.21 | 0.55 | 0.0040 | 0.0369 | 1.00 | 1.00 | 21.2 |
| human mtDNA | 16569 | 10 | 2.98 | 0.30 | 0.17 | 0.72 | 0.0131 | 0.0535 | 1.00 | 0.95 | 24.0 |

## Key numbers (main text, m = 100)

Myers, `filter`, `filter_fa`: median [10th, 90th percentile] over 100 reads, ms per query. Other columns: mean over 20 reads (Sellers 1, WFA2 5); recall columns are heuristic read recall (20 reads).

| k | k/m | Sellers | Banded | Myers | Edlib (best/read) | WFA2 (best/read) | filter | filter_fa | one-hit | two-hit | fallbacks (of 100) filter | filter_fa |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.01 | 909 | 76.4 | 70.8 [67.6, 83.5] | 44.2 | 66.8 | 0.007 [0.002, 0.018] | 0.029 [0.026, 0.034] | 0.060 (rec 1.00) | 0.036 (rec 1.00) | 0 | 0 |
| 2 | 0.02 | 943 | 93.7 | 73.0 [67.9, 87.4] | 44.0 | 91.5 | 0.017 [0.009, 0.035] | 0.034 [0.031, 0.047] | 0.077 (rec 1.00) | 0.032 (rec 1.00) | 0 | 0 |
| 3 | 0.03 | 936 | 120.3 | 70.8 [67.8, 86.6] | 39.6 | 121.7 | 0.023 [0.012, 0.040] | 0.034 [0.031, 0.048] | 0.061 (rec 1.00) | 0.032 (rec 1.00) | 0 | 0 |
| 4 | 0.04 | 933 | 133.0 | 71.4 [67.8, 82.0] | 41.5 | 156.2 | 0.021 [0.012, 0.041] | 0.036 [0.033, 0.043] | 0.063 (rec 1.00) | 0.030 (rec 1.00) | 0 | 0 |
| 6 | 0.06 | 944 | 176.8 | 69.7 [65.7, 80.2] | 43.6 | 173.7 | 0.029 [0.017, 0.052] | 0.043 [0.039, 0.053] | 0.079 (rec 1.00) | 0.028 (rec 1.00) | 0 | 0 |
| 8 | 0.08 | 845 | 207.5 | 71.0 [67.8, 77.4] | 44.9 | 253.9 | 0.040 [0.023, 0.126] | 0.084 [0.067, 0.114] | 0.073 (rec 1.00) | 0.030 (rec 0.95) | 0 | 0 |
| 10 | 0.10 | 854 | 241.8 | 71.3 [67.7, 88.6] | 37.0 | 313.1 | 0.631 [0.412, 2.309] | 0.455 [0.271, 0.693] | 0.080 (rec 1.00) | 0.037 (rec 0.95) | 0 | 0 |
| 11 | 0.11 | 862 | 250.4 | 72.4 [69.0, 87.5] | 40.9 | 333.2 | 2.581 [2.039, 3.357] | 1.494 [0.815, 2.264] | 0.051 (rec 1.00) | 0.023 (rec 0.70) | 0 | 0 |
| 12 | 0.12 | 831 | 272.5 | 70.2 [66.9, 77.6] | 40.6 | 339.3 | 9.117 [3.474, 11.778] | 3.238 [2.147, 4.967] | 0.052 (rec 1.00) | 0.026 (rec 0.75) | 0 | 0 |
| 13 | 0.13 | 931 | 295.5 | 70.9 [68.2, 89.3] | 36.7 | 321.7 | 10.709 [8.711, 29.728] | 7.256 [4.659, 11.353] | 0.065 (rec 1.00) | 0.025 (rec 0.75) | 0 | 0 |
| 14 | 0.14 | 875 | 314.2 | 69.3 [66.2, 75.5] | 49.7 | 452.7 | 32.443 [22.573, 39.132] | 13.982 [10.537, 20.034] | 0.079 (rec 1.00) | 0.032 (rec 0.70) | 0 | 0 |
| 15 | 0.15 | 886 | 337.3 | 69.6 [66.5, 74.1] | 41.9 | 382.9 | 35.638 [31.004, 43.287] | 23.340 [17.584, 31.558] | 0.048 (rec 1.00) | 0.026 (rec 0.60) | 4 | 0 |
| 17 | 0.17 | 853 | 364.8 | 68.7 [65.6, 77.5] | 43.7 | 493.0 | 70.546 [67.629, 72.557] | 53.992 [40.952, 62.860] | 0.044 (rec 0.85) | 0.021 (rec 0.45) | 87 | 1 |
| 20 | 0.20 | 852 | 413.4 | 68.9 [66.5, 71.9] | 38.3 | 564.3 | 71.311 [69.181, 78.833] | 72.574 [69.321, 77.082] | 0.052 (rec 0.95) | 0.022 (rec 0.40) | 100 | 99 |
| 30 | 0.30 | 887 | 593.4 | 70.3 [66.8, 77.7] | 55.5 | 766.1 | 70.945 [68.392, 72.966] | 71.556 [68.185, 82.860] | 0.044 (rec 0.30) | 0.023 (rec 0.00) | 100 | 100 |

## Instructions per cell and branch behaviour: banded vs Myers

Hardware performance counters: perf installed = False, CPU PMU device present = False. In this VM neither is true, so hardware counters could not be used. The numbers below come from `valgrind --tool=cachegrind --branch-sim=yes`: instruction counts are exact executed-instruction counts; branch mispredictions come from cachegrind's simulated predictor and are indicative only. They do not show cache misses or real timing.

Text: first 150,000 bp of E. coli 536; m = 100; per pattern (mean of 5). Myers cells = m * n DP cells represented (64 per word operation); banded cells = cells actually computed.

| k | algorithm | instructions per text char | instructions per cell | cells per text char | branches per text char | simulated mispredicts per 1000 chars | simulated mispredict rate |
|---|---|---|---|---|---|---|---|
| 5 | banded | 229.8 | 21.43 | 10.7 | 16.29 | 2024.7 | 12.43% |
| 5 | myers | 98.7 | 0.99 | 100.0 | 9.01 | 621.8 | 6.90% |
| 20 | banded | 754.3 | 18.92 | 39.9 | 45.42 | 1952.4 | 4.30% |
| 20 | myers | 98.8 | 0.99 | 100.0 | 9.01 | 639.3 | 7.10% |

At k = 5: banded executes 2.3x the instructions of Myers per text character while computing 9x fewer cells; instructions per cell are 21.4 (banded) against 0.99 (Myers).

At k = 20: banded executes 7.6x the instructions of Myers per text character while computing 3x fewer cells; instructions per cell are 18.9 (banded) against 0.99 (Myers).


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

Predicted crossover with a Myers scan (cost-model constants, formula hits): uniform-text formula k = 17.0 (k/m = 0.170); spectrum-corrected formula k = 17.0 (k/m = 0.170).

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

Predicted crossover with a Myers scan (cost-model constants, formula hits): uniform-text formula k = 17.0 (k/m = 0.170); spectrum-corrected formula k = 17.0 (k/m = 0.170).

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

Predicted crossover with a Myers scan (cost-model constants, formula hits): uniform-text formula k = 17.0 (k/m = 0.170); spectrum-corrected formula k = 16.7 (k/m = 0.167).

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
| 10 | 39 | 0.2652 | 0.1131 | 0.9199 | 0.0064 | 0.0431 | 5 | 2 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 110 | 0.2327 | 0.1187 | 0.9394 | 0.1536 | 0.1576 | 135 | 73 | 0 / 200 | 1.00 | 0.96 | 1.00 | 0.98 |
| 20 | 177 | 0.2290 | 0.2288 | 0.9626 | 0.3539 | 0.3775 | 2082 | 680 | 200 / 200 | 0.99 | 0.85 | 1.00 | 0.93 |
| 25 | 190 | 0.2719 | 0.1491 | 0.9337 | 0.8315 | 0.4996 | 9006 | 2786 | 200 / 200 | 0.99 | 0.82 | 1.00 | 0.89 |
| 30 | 195 | 0.2295 | 0.1562 | 0.9617 | 1.0568 | 0.8340 | 10743 | 7682 | 200 / 200 | 0.98 | 0.80 | 0.99 | 0.86 |
| 40 | 200 | 0.2320 | 0.1538 | 0.9205 | 4.4296 | 2.6464 | 51031 | 28845 | 200 / 200 | 0.96 | 0.79 | 0.77 | 0.65 |

### mt_orang segments -> mt_human text, m = 200 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 29 (14.5% of m), 10th to 90th percentile 20 to 41 (10.0% to 20.5%), min 4, max 93; 0 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 3 | 0.4600 | 0.1174 | 1.7791 | 0.0073 | 0.0726 | 2 | 0 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 11 | 0.5815 | 0.1531 | 1.9333 | 0.0089 | 0.1328 | 4 | 1 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 20 | 22 | 0.4558 | 0.1143 | 1.7449 | 0.0370 | 0.1168 | 10 | 3 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 25 | 60 | 0.4532 | 0.1320 | 1.7155 | 0.2968 | 0.2069 | 72 | 16 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 30 | 129 | 0.4757 | 0.2266 | 1.6163 | 0.4631 | 0.4721 | 265 | 105 | 66 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 40 | 179 | 0.4485 | 0.1655 | 1.6421 | 0.6913 | 0.8433 | 4073 | 1133 | 200 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |

### mt_human segments -> mt_orang text, m = 100 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 15 (15.0% of m), 10th to 90th percentile 8 to 24 (8.0% to 24.0%), min 0, max 40; 1 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 43 | 0.2251 | 0.1118 | 1.0129 | 0.0073 | 0.0353 | 5 | 2 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 111 | 0.2134 | 0.1157 | 0.9118 | 0.1518 | 0.1353 | 136 | 69 | 0 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |
| 20 | 163 | 0.2282 | 0.1267 | 0.9553 | 0.3517 | 0.3443 | 2067 | 664 | 200 / 200 | 1.00 | 0.86 | 1.00 | 0.94 |
| 25 | 182 | 0.2323 | 0.1271 | 0.9568 | 0.9892 | 0.4751 | 8979 | 2746 | 200 / 200 | 0.99 | 0.79 | 1.00 | 0.89 |
| 30 | 188 | 0.2300 | 0.1509 | 1.0040 | 1.0015 | 0.8719 | 10736 | 7541 | 200 / 200 | 0.99 | 0.78 | 1.00 | 0.85 |
| 40 | 200 | 0.2313 | 0.1997 | 1.0174 | 4.5348 | 2.6540 | 50822 | 28455 | 200 / 200 | 0.95 | 0.73 | 0.77 | 0.63 |

### mt_human segments -> mt_orang text, m = 200 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 29 (14.5% of m), 10th to 90th percentile 21 to 41 (10.5% to 20.5%), min 6, max 85; 0 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 3 | 0.4744 | 0.1375 | 1.9522 | 0.0071 | 0.0695 | 2 | 0 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 11 | 0.4495 | 0.1301 | 1.7999 | 0.0086 | 0.0879 | 4 | 1 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 20 | 19 | 0.4875 | 0.1295 | 1.8229 | 0.0366 | 0.1193 | 10 | 3 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 25 | 47 | 0.4607 | 0.2319 | 2.0027 | 0.2724 | 0.1882 | 72 | 15 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 30 | 115 | 0.5078 | 0.1763 | 1.8609 | 0.4863 | 0.5308 | 264 | 99 | 73 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 40 | 178 | 0.4897 | 0.1737 | 1.8829 | 0.7285 | 0.7168 | 4079 | 1099 | 200 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |

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

Pure filters (fallback off), mean over reads; windows = merged windows verified; lambda = (c_win + c_bp B (m + 2k)) / c_hit with the fitted constants (c_hit = 109.7, c_win = 219, c_bp = 7.28 ns), e.g. lambda = 17.9 at m = 100, k = 10. All methods return the same matches as Myers (asserted).

### Stress scenarios (E. coli 536 + about 1% microsatellite/homopolymer arrays, m = 100, 20 reads)

| scenario | k | algorithm | hits | windows | bp verified | query ms | seed selection ms |
|---|---|---|---|---|---|---|---|
| uniform | 6 | equal split | 17 | 14 | 1590 | 0.028 | 0.000 |
| uniform | 6 | fa: min hits | 4 | 2 | 241 | 0.049 | 0.042 |
| uniform | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 4 | 2 | 241 | 0.136 | 0.130 |
| uniform | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 4 | 2 | 241 | 0.160 | 0.153 |
| uniform | 10 | equal split | 463 | 448 | 54540 | 0.872 | 0.000 |
| uniform | 10 | fa: min hits | 186 | 181 | 22008 | 0.452 | 0.059 |
| uniform | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 186 | 181 | 22008 | 0.528 | 0.162 |
| uniform | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 186 | 181 | 22008 | 0.487 | 0.137 |
| uniform | 15 | equal split | 25050 | 11833 | 2165625 | 33.630 | 0.000 |
| uniform | 15 | fa: min hits | 15693 | 9462 | 1522426 | 23.500 | 0.083 |
| uniform | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16012 | 9568 | 1548422 | 26.420 | 0.668 |
| uniform | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 15693 | 9462 | 1522426 | 23.081 | 0.174 |
| clustered | 6 | equal split | 22 | 16 | 1793 | 0.034 | 0.000 |
| clustered | 6 | fa: min hits | 8 | 3 | 384 | 0.051 | 0.042 |
| clustered | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 8 | 3 | 384 | 0.160 | 0.147 |
| clustered | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 8 | 3 | 384 | 0.153 | 0.144 |
| clustered | 10 | equal split | 3256 | 513 | 62756 | 1.101 | 0.000 |
| clustered | 10 | fa: min hits | 2446 | 376 | 47397 | 0.884 | 0.050 |
| clustered | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 2616 | 228 | 28192 | 3.204 | 2.640 |
| clustered | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 2446 | 376 | 47392 | 1.042 | 0.153 |
| clustered | 15 | equal split | 28859 | 11635 | 2235284 | 36.641 | 0.000 |
| clustered | 15 | fa: min hits | 18291 | 9384 | 1530393 | 25.585 | 0.086 |
| clustered | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 19064 | 9771 | 1618158 | 31.182 | 2.581 |
| clustered | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 18291 | 9384 | 1530393 | 24.866 | 0.198 |
| indel_heavy | 6 | equal split | 17 | 14 | 1528 | 0.029 | 0.000 |
| indel_heavy | 6 | fa: min hits | 5 | 2 | 275 | 0.049 | 0.042 |
| indel_heavy | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 5 | 2 | 264 | 0.129 | 0.122 |
| indel_heavy | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 5 | 2 | 264 | 0.133 | 0.126 |
| indel_heavy | 10 | equal split | 679 | 655 | 78334 | 1.275 | 0.000 |
| indel_heavy | 10 | fa: min hits | 248 | 242 | 28661 | 0.488 | 0.048 |
| indel_heavy | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 248 | 242 | 28661 | 0.619 | 0.139 |
| indel_heavy | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 248 | 242 | 28661 | 0.585 | 0.129 |
| indel_heavy | 15 | equal split | 33642 | 11860 | 2568773 | 42.434 | 0.000 |
| indel_heavy | 15 | fa: min hits | 16342 | 10085 | 1668047 | 25.836 | 0.081 |
| indel_heavy | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16342 | 10085 | 1668047 | 27.737 | 0.216 |
| indel_heavy | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 16342 | 10085 | 1668047 | 27.896 | 0.209 |
| lowcomplex_edge | 6 | equal split | 6878 | 21 | 12819 | 0.640 | 0.000 |
| lowcomplex_edge | 6 | fa: min hits | 12 | 10 | 1148 | 0.060 | 0.040 |
| lowcomplex_edge | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 12 | 10 | 1148 | 11.767 | 11.744 |
| lowcomplex_edge | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 12 | 10 | 1148 | 0.314 | 0.287 |
| lowcomplex_edge | 10 | equal split | 12693 | 531 | 75968 | 1.707 | 0.000 |
| lowcomplex_edge | 10 | fa: min hits | 2839 | 1536 | 193687 | 3.242 | 0.054 |
| lowcomplex_edge | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 5662 | 295 | 42867 | 8.767 | 7.833 |
| lowcomplex_edge | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 2859 | 1515 | 191528 | 3.635 | 0.192 |
| lowcomplex_edge | 15 | equal split | 38737 | 10889 | 1938189 | 30.480 | 0.000 |
| lowcomplex_edge | 15 | fa: min hits | 26362 | 9964 | 1696660 | 26.833 | 0.073 |
| lowcomplex_edge | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 28203 | 10940 | 1934724 | 36.549 | 4.238 |
| lowcomplex_edge | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 26362 | 9964 | 1696660 | 28.156 | 0.259 |
| lowcomplex_inside | 6 | equal split | 17734 | 14 | 11878 | 0.879 | 0.000 |
| lowcomplex_inside | 6 | fa: min hits | 7111 | 140 | 26020 | 0.725 | 0.031 |
| lowcomplex_inside | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 7721 | 35 | 13841 | 23.443 | 22.776 |
| lowcomplex_inside | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 7133 | 121 | 23833 | 1.258 | 0.507 |
| lowcomplex_inside | 10 | equal split | 18046 | 535 | 85941 | 2.102 | 0.000 |
| lowcomplex_inside | 10 | fa: min hits | 11129 | 501 | 71955 | 1.603 | 0.042 |
| lowcomplex_inside | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 12009 | 326 | 50918 | 14.753 | 13.238 |
| lowcomplex_inside | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 11150 | 497 | 71543 | 2.291 | 0.603 |
| lowcomplex_inside | 15 | equal split | 51356 | 7508 | 1461597 | 25.737 | 0.000 |
| lowcomplex_inside | 15 | fa: min hits | 33393 | 6438 | 1121207 | 19.895 | 0.066 |
| lowcomplex_inside | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 36042 | 7451 | 1354083 | 30.815 | 7.343 |
| lowcomplex_inside | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 33454 | 6506 | 1130474 | 20.759 | 0.396 |

### Unmodified E. coli 536, uniform edits, m = 100, 100 reads (regression check)

| k | algorithm | hits | windows | bp verified | query ms | seed selection ms |
|---|---|---|---|---|---|---|
| 4 | equal split | 16 | 13 | 1478 | 0.035 | 0.000 |
| 4 | fa: min hits | 3 | 2 | 177 | 0.049 | 0.043 |
| 4 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 3 | 2 | 177 | 0.172 | 0.166 |
| 4 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 3 | 2 | 177 | 0.153 | 0.149 |
| 8 | equal split | 34 | 30 | 3479 | 0.061 | 0.000 |
| 8 | fa: min hits | 16 | 13 | 1509 | 0.071 | 0.043 |
| 8 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16 | 13 | 1504 | 0.205 | 0.169 |
| 8 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 16 | 13 | 1504 | 0.191 | 0.156 |
| 10 | equal split | 493 | 479 | 57545 | 0.917 | 0.000 |
| 10 | fa: min hits | 215 | 210 | 25137 | 0.453 | 0.049 |
| 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 215 | 210 | 25136 | 0.617 | 0.159 |
| 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 215 | 210 | 25136 | 0.609 | 0.150 |
| 12 | equal split | 4654 | 4050 | 538390 | 8.589 | 0.000 |
| 12 | fa: min hits | 1760 | 1665 | 211379 | 3.594 | 0.068 |
| 12 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 1760 | 1665 | 211379 | 3.686 | 0.178 |
| 12 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 1760 | 1665 | 211379 | 3.943 | 0.187 |
| 15 | equal split | 28807 | 12005 | 2422582 | 38.869 | 0.000 |
| 15 | fa: min hits | 14685 | 9594 | 1551347 | 25.044 | 0.087 |
| 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 14685 | 9594 | 1551347 | 24.462 | 0.193 |
| 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 14685 | 9594 | 1551347 | 25.212 | 0.199 |

**Before/after on the case that verified more bases (lowcomplex_edge, k = 10):** equal split 12693 hits, 531 windows, 75968 bp, 1.71 ms; min-hits DP 2839 hits, 1536 windows, 193687 bp, 3.24 ms; hits + lambda * windows DP 5662 hits, 295 windows, 42867 bp, 8.77 ms (seed selection 7.83 ms); same with hit-list cap 1000: 2859 hits, 1515 windows, 191528 bp, 3.64 ms (selection 0.19 ms).
`filter_fa_w` versus `filter_fa` over all 20 settings: more than 5% faster in 0, more than 5% slower in 18, within 5% in 2.
`filter_fa_w1k` versus `filter_fa` over all 20 settings: more than 5% faster in 0, more than 5% slower in 15, within 5% in 5.


## Peak memory per method (measured)

E. coli 536 (n = 4,938,920), m = 100, k = 10, 5 reads; peak resident set size (VmHWM). working memory = peak RSS minus the steady RSS after the text and reads are loaded and the heap trimmed (8.4 MB); methods whose peak does not exceed the transient peak of loading the file (about 20 MB, the file is read line by line into a string) are shown as negligible. Edlib and WFA2 run in a Python process, so their baseline differs.

| method | peak RSS (MB) | working memory of the method (MB) | expected |
|---|---|---|---|
| noop | 20.2 | negligible (below the loading peak) | text only (n bytes) plus read buffers |
| sellers | 20.2 | negligible (below the loading peak) | O(m) |
| banded | 20.2 | negligible (below the loading peak) | O(m) |
| myers | 20.2 | negligible (below the loading peak) | O(B sigma) |
| myers_pruned | 20.1 | negligible (below the loading peak) | O(B sigma) |
| lv | 50.8 | 42.4 | 2 arrays of n + 2k ints (8n bytes) |
| filter | 63.9 | 55.5 | index 4(4^11+1) + 4n bytes (37 MB) plus build temporaries (counter copy 17 MB, flags 5 MB) |
| filter_fa | 63.9 | 55.5 | same index |
| filter_fa_w | 63.9 | 55.5 | same index |
| blast1 | 63.9 | 55.5 | same index |
| blast2 | 63.9 | 55.5 | same index |
| map | 63.9 | 55.5 | same index + O(m (m + k)) per cluster |
| edlib (Python process) | 39.3 | negligible (below the loading peak) | library; measured as ru_maxrss increase over loading the text |

