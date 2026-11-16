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
