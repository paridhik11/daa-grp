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
