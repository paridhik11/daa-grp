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
| stress low-complexity edge | 6 | 0.477 (0.000) | 0.042 (0.005) | 11.558 (11.520) | 0.680 (0.646) | 0.700 (0.664) | 11.349 (11.310) | 0.042 (0.005) |
| stress low-complexity edge | 10 | 1.503 (0.000) | 2.603 (0.015) | 9.560 (8.643) | 1.495 (0.549) | 1.517 (0.547) | 9.634 (8.680) | 1.435 (0.014) |
| stress low-complexity edge | 15 | 34.572 (0.000) | 25.893 (0.025) | 33.373 (4.513) | 29.097 (0.361) | 28.239 (0.027) | 27.246 (0.026) | 25.000 (0.023) |
| stress low-complexity inside | 6 | 0.771 (0.000) | 0.569 (0.006) | 24.582 (23.953) | 1.832 (1.288) | 1.904 (1.335) | 23.506 (22.909) | 0.717 (0.007) |
| stress low-complexity inside | 10 | 1.889 (0.000) | 1.580 (0.011) | 17.498 (16.085) | 2.288 (0.937) | 2.264 (0.911) | 16.611 (15.272) | 1.482 (0.010) |
| stress low-complexity inside | 15 | 26.560 (0.000) | 20.691 (0.020) | 36.679 (10.218) | 25.715 (0.627) | 19.285 (0.018) | 19.476 (0.018) | 20.019 (0.018) |
| E. coli uniform | 6 | 0.033 (0.000) | 0.013 (0.006) | 0.132 (0.124) | 0.053 (0.046) | 0.013 (0.006) | 0.014 (0.006) | 0.012 (0.005) |
| E. coli uniform | 10 | 0.928 (0.000) | 0.425 (0.012) | 0.725 (0.182) | 0.500 (0.071) | 0.428 (0.013) | 0.427 (0.013) | 0.437 (0.014) |
| E. coli uniform | 15 | 39.063 (0.000) | 24.400 (0.026) | 25.075 (0.191) | 24.650 (0.119) | 25.045 (0.027) | 25.397 (0.027) | 25.329 (0.026) |
| stress clustered | 6 | 0.046 (0.000) | 0.023 (0.006) | 0.365 (0.348) | 0.080 (0.064) | 0.037 (0.020) | 0.235 (0.218) | 0.024 (0.007) |
| stress clustered | 10 | 1.108 (0.000) | 0.510 (0.013) | 1.251 (0.740) | 0.634 (0.105) | 0.534 (0.039) | 1.054 (0.559) | 0.536 (0.014) |
| stress clustered | 15 | 38.551 (0.000) | 25.297 (0.028) | 26.564 (0.852) | 26.613 (0.167) | 25.637 (0.026) | 26.377 (0.029) | 26.775 (0.027) |
