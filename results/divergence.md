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
