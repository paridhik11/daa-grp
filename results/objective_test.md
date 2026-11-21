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
