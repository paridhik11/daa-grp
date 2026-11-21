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
