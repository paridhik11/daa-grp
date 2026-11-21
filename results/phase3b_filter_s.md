## Phase 3B, item 4: generalised pigeonhole (s required hits in a diagonal band of width k)

k + s disjoint seeds; a candidate needs s distinct seeds hit within a band of k diagonals; s = 1 is the classic filter. Pure filters (cost-model fallback off), mean per read; `windows` = merged candidate windows verified, `bp` = bases verified, `gen` = time to build candidates (sort and band pass). Every filter's matches are asserted equal to Myers on every read.


### E. coli 536, m = 100 (10% error), k = 10

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 75.678 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 1.238 | 508 | 493 | 59307 | 0.029 | 0.000 |
| filter_fa (fa, s=1) | 0.588 | 220 | 215 | 25771 | 0.015 | 0.019 |
| filter_s2 (equal, s=2) | 0.082 | 1368 | 5 | 678 | 0.069 | 0.000 |
| filter_s3 (equal, s=3) | 0.274 | 4813 | 1 | 145 | 0.268 | 0.000 |
| filter_fs2 (fa, s=2) | 0.049 | 678 | 2 | 241 | 0.031 | 0.011 |
| filter_fs3 (fa, s=3) | 0.103 | 1735 | 1 | 130 | 0.089 | 0.011 |

### E. coli 536, m = 100 (12% error), k = 12

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 75.303 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 8.434 | 4621 | 4023 | 532929 | 0.298 | 0.000 |
| filter_fa (fa, s=1) | 3.631 | 1802 | 1702 | 215592 | 0.109 | 0.018 |
| filter_s2 (equal, s=2) | 0.984 | 8147 | 231 | 29942 | 0.528 | 0.000 |
| filter_s3 (equal, s=3) | 1.439 | 21128 | 27 | 3548 | 1.386 | 0.000 |
| filter_fs2 (fa, s=2) | 0.329 | 3986 | 40 | 5187 | 0.224 | 0.013 |
| filter_fs3 (fa, s=3) | 0.523 | 8188 | 4 | 467 | 0.527 | 0.015 |

### E. coli 536, m = 100 (15% error), k = 15

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 75.254 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 42.950 | 33282 | 11915 | 2581677 | 2.818 | 0.000 |
| filter_fa (fa, s=1) | 25.572 | 14980 | 9720 | 1573613 | 1.168 | 0.029 |
| filter_s2 (equal, s=2) | 25.081 | 66872 | 7123 | 1266139 | 5.625 | 0.000 |
| filter_s3 (equal, s=3) | 13.615 | 98311 | 2474 | 368362 | 8.395 | 0.000 |
| filter_fs2 (fa, s=2) | 5.981 | 26365 | 1839 | 262385 | 2.040 | 0.027 |
| filter_fs3 (fa, s=3) | 3.972 | 44034 | 339 | 47615 | 3.288 | 0.021 |

### E. coli 536, m = 100 (17% error), k = 17

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 79.607 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 77.715 | 97109 | 7007 | 4550851 | 9.322 | 0.000 |
| filter_fa (fa, s=1) | 54.495 | 43480 | 12610 | 3303363 | 3.775 | 0.032 |
| filter_s2 (equal, s=2) | 48.431 | 120521 | 10765 | 2578861 | 9.936 | 0.000 |
| filter_s3 (equal, s=3) | 52.271 | 232017 | 4772 | 2163209 | 20.501 | 0.000 |
| filter_fs2 (fa, s=2) | 27.947 | 67827 | 7723 | 1384864 | 5.573 | 0.030 |
| filter_fs3 (fa, s=3) | 17.285 | 102493 | 3380 | 549434 | 8.150 | 0.028 |

### E. coli 536, m = 100 (20% error), k = 20

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 78.797 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 114.097 | 400572 | 521 | 4913075 | 42.230 | 0.000 |
| filter_fa (fa, s=1) | 85.346 | 150806 | 2859 | 4813202 | 13.523 | 0.035 |
| filter_s2 (equal, s=2) | 115.084 | 463463 | 222 | 4932888 | 43.233 | 0.000 |
| filter_s3 (equal, s=3) | 119.190 | 484897 | 2199 | 4837082 | 48.602 | 0.000 |
| filter_fs2 (fa, s=2) | 89.646 | 213659 | 6764 | 4387587 | 19.199 | 0.035 |
| filter_fs3 (fa, s=3) | 83.850 | 292128 | 9138 | 3761239 | 27.406 | 0.034 |

### E. coli 536, m = 100 (25% error), k = 25

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 73.376 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 261.709 | 1605638 | 1 | 4938919 | 180.567 | 0.000 |
| filter_fa (fa, s=1) | 141.522 | 633545 | 2 | 4938912 | 70.401 | 0.041 |
| filter_s2 (equal, s=2) | 302.836 | 2089283 | 2 | 4938903 | 237.039 | 0.000 |
| filter_s3 (equal, s=3) | 308.442 | 2267775 | 1 | 4938920 | 237.262 | 0.000 |
| filter_fs2 (fa, s=2) | 154.237 | 809200 | 3 | 4938880 | 81.422 | 0.038 |
| filter_fs3 (fa, s=3) | 180.331 | 1014714 | 6 | 4938806 | 106.887 | 0.040 |

### mt_orang segments vs mt_human (real homologs), m = 100, k = 15

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.249 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.160 | 133 | 45 | 10074 | 0.006 | 0.000 |
| filter_fa (fa, s=1) | 0.116 | 71 | 38 | 6528 | 0.006 | 0.017 |
| filter_s2 (equal, s=2) | 0.176 | 474 | 40 | 10018 | 0.024 | 0.000 |
| filter_s3 (equal, s=3) | 0.078 | 503 | 20 | 3345 | 0.025 | 0.000 |
| filter_fs2 (fa, s=2) | 0.043 | 123 | 10 | 1496 | 0.006 | 0.013 |
| filter_fs3 (fa, s=3) | 0.030 | 202 | 4 | 492 | 0.008 | 0.012 |

### mt_orang segments vs mt_human (real homologs), m = 100, k = 20

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.246 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.360 | 2046 | 1 | 16569 | 0.115 | 0.000 |
| filter_fa (fa, s=1) | 0.353 | 665 | 5 | 16429 | 0.033 | 0.014 |
| filter_s2 (equal, s=2) | 0.467 | 2139 | 1 | 16562 | 0.138 | 0.000 |
| filter_s3 (equal, s=3) | 0.486 | 2239 | 3 | 16492 | 0.166 | 0.000 |
| filter_fs2 (fa, s=2) | 0.321 | 931 | 13 | 15854 | 0.052 | 0.016 |
| filter_fs3 (fa, s=3) | 0.328 | 1273 | 21 | 14905 | 0.078 | 0.016 |

### mt_orang segments vs mt_human (real homologs), m = 100, k = 25

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.253 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 1.082 | 8985 | 1 | 16569 | 0.658 | 0.000 |
| filter_fa (fa, s=1) | 0.573 | 2748 | 1 | 16569 | 0.183 | 0.023 |
| filter_s2 (equal, s=2) | 0.980 | 9344 | 1 | 16569 | 0.847 | 0.000 |
| filter_s3 (equal, s=3) | 1.200 | 9680 | 1 | 16569 | 0.896 | 0.000 |
| filter_fs2 (fa, s=2) | 0.495 | 3438 | 1 | 16569 | 0.231 | 0.022 |
| filter_fs3 (fa, s=3) | 0.606 | 4282 | 1 | 16569 | 0.288 | 0.020 |

### mt_orang segments vs mt_human (real homologs), m = 100, k = 30

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.283 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 1.009 | 10716 | 1 | 16569 | 0.987 | 0.000 |
| filter_fa (fa, s=1) | 1.005 | 7622 | 1 | 16569 | 0.550 | 0.027 |
| filter_s2 (equal, s=2) | 1.145 | 11050 | 1 | 16569 | 0.983 | 0.000 |
| filter_s3 (equal, s=3) | 1.043 | 11385 | 1 | 16569 | 0.792 | 0.000 |
| filter_fs2 (fa, s=2) | 0.882 | 8935 | 1 | 16569 | 0.808 | 0.038 |
| filter_fs3 (fa, s=3) | 1.007 | 10351 | 1 | 16569 | 0.716 | 0.021 |

### mt_orang segments vs mt_human (real homologs), m = 200, k = 30

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.602 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.510 | 267 | 7 | 16060 | 0.016 | 0.000 |
| filter_fa (fa, s=1) | 0.543 | 106 | 21 | 12394 | 0.006 | 0.037 |
| filter_s2 (equal, s=2) | 0.362 | 276 | 20 | 11305 | 0.019 | 0.000 |
| filter_s3 (equal, s=3) | 0.155 | 284 | 13 | 4392 | 0.014 | 0.000 |
| filter_fs2 (fa, s=2) | 0.224 | 145 | 15 | 5452 | 0.007 | 0.035 |
| filter_fs3 (fa, s=3) | 0.147 | 194 | 6 | 1970 | 0.012 | 0.058 |

### mt_orang segments vs mt_human (real homologs), m = 200, k = 40

| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |
|---|---|---|---|---|---|---|
| Myers scan | 0.492 | 0 | 0 | 0 | 0.000 | 0.000 |
| filter (equal, s=1) | 0.750 | 4094 | 1 | 16569 | 0.252 | 0.000 |
| filter_fa (fa, s=1) | 0.607 | 1150 | 1 | 16569 | 0.065 | 0.048 |
| filter_s2 (equal, s=2) | 0.840 | 4181 | 1 | 16569 | 0.298 | 0.000 |
| filter_s3 (equal, s=3) | 0.813 | 4284 | 1 | 16569 | 0.354 | 0.000 |
| filter_fs2 (fa, s=2) | 0.631 | 1369 | 1 | 16568 | 0.080 | 0.048 |
| filter_fs3 (fa, s=3) | 0.650 | 1620 | 1 | 16566 | 0.098 | 0.050 |
