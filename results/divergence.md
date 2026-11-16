## Real divergence: mt_orang segments against mt_human, and the reverse

Both are real homologous mitochondrial genomes (human, 16,569 bp; orangutan, 16,499 bp), so the edits are real evolutionary differences, not simulated ones. Segments are random substrings of the query genome; the text is the other genome (forward strand). The true edit distance is the minimum over all text positions, computed with Myers and verified with Sellers.


### mt_orang segments -> mt_human text, m = 100 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 15 (15.0% of m), 10th to 90th percentile 8 to 22 (8.0% to 22.0%), min 1, max 39; 0 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 39 | 0.2520 | 0.1141 | 0.9031 | 0.0070 | 0.0261 | 5 | 2 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 110 | 0.2362 | 0.1352 | 1.1333 | 0.1600 | 0.1452 | 135 | 73 | 1 / 200 | 1.00 | 0.96 | 1.00 | 0.98 |
| 20 | 177 | 0.2622 | 0.1903 | 0.9258 | 0.4744 | 0.3794 | 2082 | 680 | 200 / 200 | 0.99 | 0.85 | 1.00 | 0.93 |
| 25 | 190 | 0.2324 | 0.1367 | 0.9269 | 0.9847 | 0.4525 | 9006 | 2786 | 200 / 200 | 0.99 | 0.82 | 1.00 | 0.89 |
| 30 | 195 | 0.2650 | 0.2398 | 0.9463 | 1.1399 | 0.8230 | 10743 | 7682 | 200 / 200 | 0.98 | 0.80 | 0.99 | 0.86 |
| 40 | 200 | 0.2272 | 0.1565 | 0.9210 | 4.4420 | 2.5693 | 51031 | 28845 | 200 / 200 | 0.96 | 0.79 | 0.77 | 0.65 |

### mt_orang segments -> mt_human text, m = 200 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 29 (14.5% of m), 10th to 90th percentile 20 to 41 (10.0% to 20.5%), min 4, max 93; 0 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 3 | 0.4661 | 0.1767 | 1.7427 | 0.0071 | 0.0572 | 2 | 0 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 11 | 0.4744 | 0.2164 | 1.8212 | 0.0094 | 0.0556 | 4 | 1 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 20 | 22 | 0.4855 | 0.1210 | 1.7703 | 0.0386 | 0.0754 | 10 | 3 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 25 | 60 | 0.4733 | 0.1533 | 1.8442 | 0.2898 | 0.1461 | 72 | 16 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 30 | 129 | 0.4754 | 0.2873 | 2.1236 | 0.4706 | 0.4635 | 265 | 105 | 168 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 40 | 179 | 0.4745 | 0.1849 | 1.9012 | 0.8672 | 0.6899 | 4073 | 1133 | 200 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |

### mt_human segments -> mt_orang text, m = 100 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 15 (15.0% of m), 10th to 90th percentile 8 to 24 (8.0% to 24.0%), min 0, max 40; 1 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 43 | 0.2291 | 0.1130 | 0.9812 | 0.0087 | 0.0272 | 5 | 2 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 111 | 0.2381 | 0.1447 | 0.9955 | 0.1628 | 0.1356 | 136 | 69 | 5 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |
| 20 | 163 | 0.3033 | 0.1395 | 1.3117 | 0.4506 | 0.3057 | 2067 | 664 | 200 / 200 | 1.00 | 0.86 | 1.00 | 0.94 |
| 25 | 182 | 0.2703 | 0.1311 | 0.9300 | 0.8812 | 0.4430 | 8979 | 2746 | 200 / 200 | 0.99 | 0.79 | 1.00 | 0.89 |
| 30 | 188 | 0.2402 | 0.1554 | 0.9318 | 0.9603 | 0.9108 | 10736 | 7541 | 200 / 200 | 0.99 | 0.78 | 1.00 | 0.85 |
| 40 | 200 | 0.2289 | 0.1566 | 1.0615 | 4.6212 | 2.5098 | 50822 | 28455 | 200 / 200 | 0.95 | 0.73 | 0.77 | 0.63 |

### mt_human segments -> mt_orang text, m = 200 (200 segments)

True best edit distance d* (Sellers-verified, 0 verification failures): median 29 (14.5% of m), 10th to 90th percentile 21 to 41 (10.5% to 20.5%), min 6, max 85; 0 segments match exactly.

| k | segments with d* <= k | Myers ms | Edlib ms | WFA2 ms | filter (pure) ms | filter_fa (pure) ms | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 3 | 0.5049 | 0.1179 | 1.8834 | 0.0103 | 0.0664 | 2 | 0 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15 | 11 | 0.4661 | 0.1169 | 1.7901 | 0.0095 | 0.0618 | 4 | 1 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 20 | 19 | 0.5087 | 0.1166 | 1.7817 | 0.0388 | 0.0809 | 10 | 3 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 25 | 47 | 0.4793 | 0.1315 | 2.1723 | 0.2863 | 0.1520 | 72 | 15 | 0 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 30 | 115 | 0.5129 | 0.2034 | 1.8192 | 0.5645 | 0.4839 | 264 | 99 | 166 / 200 | 1.00 | 1.00 | 1.00 | 1.00 |
| 40 | 178 | 0.4772 | 0.1821 | 1.8619 | 0.7364 | 0.7561 | 4079 | 1099 | 200 / 200 | 1.00 | 0.99 | 1.00 | 1.00 |

Times are mean ms per segment on a 16.5 kbp text (all exact methods agree with Sellers on every segment and k; Edlib and WFA2 return the best match per read and include Python overhead). With a text this small the index costs about 20 ms to build (clearing the 4^11 table), so the filter is only worthwhile over many reads. "pure" = fallback disabled; the fallback column counts how many segments the cost-model rule sent to a full Myers scan.
