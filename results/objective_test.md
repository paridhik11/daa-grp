## Seed-selection objective: hits versus hits + lambda * windows

Pure filters (fallback off), mean over reads; windows = merged windows verified; lambda = (c_win + c_bp B (m + 2k)) / c_hit with the fitted constants (c_hit = 82.1, c_win = 404, c_bp = 5.10 ns), e.g. lambda = 19.8 at m = 100, k = 10. All methods return the same matches as Myers (asserted).

### Stress scenarios (E. coli 536 + about 1% microsatellite/homopolymer arrays, m = 100, 20 reads)

| scenario | k | algorithm | hits | windows | bp verified | query ms | seed selection ms |
|---|---|---|---|---|---|---|---|
| uniform | 6 | equal split | 17 | 14 | 1590 | 0.030 | 0.000 |
| uniform | 6 | fa: min hits | 4 | 2 | 241 | 0.048 | 0.042 |
| uniform | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 4 | 2 | 241 | 0.167 | 0.158 |
| uniform | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 4 | 2 | 241 | 0.141 | 0.134 |
| uniform | 10 | equal split | 463 | 448 | 54540 | 0.879 | 0.000 |
| uniform | 10 | fa: min hits | 186 | 181 | 22008 | 0.499 | 0.073 |
| uniform | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 186 | 181 | 22008 | 0.503 | 0.141 |
| uniform | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 186 | 181 | 22008 | 0.489 | 0.141 |
| uniform | 15 | equal split | 25050 | 11833 | 2165625 | 38.922 | 0.000 |
| uniform | 15 | fa: min hits | 15693 | 9462 | 1522426 | 25.902 | 0.089 |
| uniform | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16012 | 9568 | 1548422 | 24.287 | 0.706 |
| uniform | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 15693 | 9462 | 1522426 | 23.817 | 0.187 |
| clustered | 6 | equal split | 22 | 16 | 1793 | 0.032 | 0.000 |
| clustered | 6 | fa: min hits | 8 | 3 | 384 | 0.054 | 0.045 |
| clustered | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 8 | 3 | 389 | 0.190 | 0.175 |
| clustered | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 8 | 3 | 389 | 0.196 | 0.180 |
| clustered | 10 | equal split | 3256 | 513 | 62756 | 1.222 | 0.000 |
| clustered | 10 | fa: min hits | 2446 | 376 | 47397 | 0.898 | 0.053 |
| clustered | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 2616 | 228 | 28192 | 3.453 | 2.876 |
| clustered | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 2446 | 376 | 47394 | 1.004 | 0.160 |
| clustered | 15 | equal split | 28859 | 11635 | 2235284 | 37.644 | 0.000 |
| clustered | 15 | fa: min hits | 18291 | 9384 | 1530393 | 28.965 | 0.096 |
| clustered | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 19064 | 9773 | 1617684 | 29.151 | 2.544 |
| clustered | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 18291 | 9384 | 1530393 | 24.334 | 0.193 |
| indel_heavy | 6 | equal split | 17 | 14 | 1528 | 0.028 | 0.000 |
| indel_heavy | 6 | fa: min hits | 5 | 2 | 275 | 0.052 | 0.043 |
| indel_heavy | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 5 | 2 | 275 | 0.129 | 0.120 |
| indel_heavy | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 5 | 2 | 275 | 0.174 | 0.163 |
| indel_heavy | 10 | equal split | 679 | 655 | 78334 | 1.623 | 0.000 |
| indel_heavy | 10 | fa: min hits | 248 | 242 | 28661 | 0.657 | 0.064 |
| indel_heavy | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 248 | 242 | 28661 | 0.779 | 0.178 |
| indel_heavy | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 248 | 242 | 28661 | 0.621 | 0.142 |
| indel_heavy | 15 | equal split | 33642 | 11860 | 2568773 | 44.964 | 0.000 |
| indel_heavy | 15 | fa: min hits | 16342 | 10085 | 1668047 | 28.830 | 0.092 |
| indel_heavy | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16342 | 10085 | 1668047 | 30.350 | 0.230 |
| indel_heavy | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 16342 | 10085 | 1668047 | 26.547 | 0.194 |
| lowcomplex_edge | 6 | equal split | 6878 | 21 | 12819 | 0.520 | 0.000 |
| lowcomplex_edge | 6 | fa: min hits | 12 | 10 | 1148 | 0.055 | 0.034 |
| lowcomplex_edge | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 12 | 10 | 1153 | 12.808 | 12.784 |
| lowcomplex_edge | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 12 | 10 | 1153 | 0.228 | 0.209 |
| lowcomplex_edge | 10 | equal split | 12693 | 531 | 75968 | 2.150 | 0.000 |
| lowcomplex_edge | 10 | fa: min hits | 2839 | 1536 | 193687 | 3.216 | 0.051 |
| lowcomplex_edge | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 5662 | 295 | 42867 | 8.765 | 7.945 |
| lowcomplex_edge | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 2859 | 1515 | 191519 | 3.300 | 0.167 |
| lowcomplex_edge | 15 | equal split | 38737 | 10889 | 1938189 | 31.724 | 0.000 |
| lowcomplex_edge | 15 | fa: min hits | 26362 | 9964 | 1696660 | 26.358 | 0.074 |
| lowcomplex_edge | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 28203 | 10940 | 1934724 | 37.314 | 4.741 |
| lowcomplex_edge | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 26362 | 9964 | 1696660 | 28.235 | 0.262 |
| lowcomplex_inside | 6 | equal split | 17734 | 14 | 11878 | 1.136 | 0.000 |
| lowcomplex_inside | 6 | fa: min hits | 7111 | 140 | 26020 | 0.727 | 0.030 |
| lowcomplex_inside | 6 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 7721 | 35 | 13848 | 24.110 | 23.545 |
| lowcomplex_inside | 6 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 7133 | 121 | 23831 | 1.157 | 0.482 |
| lowcomplex_inside | 10 | equal split | 18046 | 535 | 85941 | 2.083 | 0.000 |
| lowcomplex_inside | 10 | fa: min hits | 11129 | 501 | 71955 | 1.572 | 0.044 |
| lowcomplex_inside | 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 12009 | 325 | 50906 | 14.042 | 12.747 |
| lowcomplex_inside | 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 11150 | 497 | 71530 | 2.342 | 0.578 |
| lowcomplex_inside | 15 | equal split | 51356 | 7508 | 1461597 | 24.299 | 0.000 |
| lowcomplex_inside | 15 | fa: min hits | 33393 | 6438 | 1121207 | 18.229 | 0.064 |
| lowcomplex_inside | 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 36042 | 7451 | 1354083 | 28.958 | 7.362 |
| lowcomplex_inside | 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 33454 | 6506 | 1130474 | 19.713 | 0.378 |

### Unmodified E. coli 536, uniform edits, m = 100, 100 reads (regression check)

| k | algorithm | hits | windows | bp verified | query ms | seed selection ms |
|---|---|---|---|---|---|---|
| 4 | equal split | 16 | 13 | 1478 | 0.025 | 0.000 |
| 4 | fa: min hits | 3 | 2 | 177 | 0.039 | 0.035 |
| 4 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 3 | 2 | 177 | 0.137 | 0.131 |
| 4 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 3 | 2 | 177 | 0.132 | 0.128 |
| 8 | equal split | 34 | 30 | 3479 | 0.058 | 0.000 |
| 8 | fa: min hits | 16 | 13 | 1509 | 0.074 | 0.046 |
| 8 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 16 | 13 | 1503 | 0.153 | 0.126 |
| 8 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 16 | 13 | 1503 | 0.169 | 0.140 |
| 10 | equal split | 493 | 479 | 57545 | 0.903 | 0.000 |
| 10 | fa: min hits | 215 | 210 | 25137 | 0.492 | 0.054 |
| 10 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 215 | 210 | 25136 | 0.537 | 0.133 |
| 10 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 215 | 210 | 25136 | 0.645 | 0.162 |
| 12 | equal split | 4654 | 4050 | 538390 | 8.466 | 0.000 |
| 12 | fa: min hits | 1760 | 1665 | 211379 | 3.409 | 0.061 |
| 12 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 1760 | 1665 | 211385 | 3.841 | 0.181 |
| 12 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 1760 | 1665 | 211385 | 3.495 | 0.164 |
| 15 | equal split | 28807 | 12005 | 2422582 | 37.772 | 0.000 |
| 15 | fa: min hits | 14685 | 9594 | 1551347 | 25.641 | 0.087 |
| 15 | fa_w: min hits + lambda*windows (hit-list cap 20000) | 14685 | 9594 | 1551347 | 26.167 | 0.204 |
| 15 | fa_w, hit-list cap 1000 (tried after seeing the first result) | 14685 | 9594 | 1551347 | 24.927 | 0.195 |

**Before/after on the case that verified more bases (lowcomplex_edge, k = 10):** equal split 12693 hits, 531 windows, 75968 bp, 2.15 ms; min-hits DP 2839 hits, 1536 windows, 193687 bp, 3.22 ms; hits + lambda * windows DP 5662 hits, 295 windows, 42867 bp, 8.76 ms (seed selection 7.94 ms); same with hit-list cap 1000: 2859 hits, 1515 windows, 191519 bp, 3.30 ms (selection 0.17 ms).
`filter_fa_w` versus `filter_fa` over all 20 settings: more than 5% faster in 1, more than 5% slower in 16, within 5% in 3.
`filter_fa_w1k` versus `filter_fa` over all 20 settings: more than 5% faster in 4, more than 5% slower in 12, within 5% in 4.
