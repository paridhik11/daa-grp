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
| stress low-complexity edge | 6 | 0.478 (0.000) | 0.042 (0.005) | 11.071 (11.034) | 0.679 (0.644) | 0.655 (0.622) | 11.974 (11.929) | 0.040 (0.011) |
| stress low-complexity edge | 10 | 1.487 (0.000) | 2.551 (0.014) | 9.538 (8.640) | 1.452 (0.529) | 1.462 (0.519) | 9.602 (8.663) | 0.444 (0.025) |
| stress low-complexity edge | 15 | 32.910 (0.000) | 25.433 (0.025) | 32.193 (4.386) | 29.158 (0.362) | 24.973 (0.022) | 24.769 (0.021) | 6.552 (0.037) |
| stress low-complexity inside | 6 | 0.781 (0.000) | 0.572 (0.005) | 23.633 (23.058) | 1.824 (1.273) | 1.848 (1.301) | 23.658 (23.070) | 0.595 (0.012) |
| stress low-complexity inside | 10 | 1.784 (0.000) | 1.565 (0.010) | 16.845 (15.525) | 2.233 (0.922) | 2.306 (0.931) | 17.071 (15.706) | 1.166 (0.024) |
| stress low-complexity inside | 15 | 25.520 (0.000) | 19.708 (0.019) | 32.911 (9.011) | 24.778 (0.592) | 19.500 (0.017) | 19.524 (0.018) | 9.654 (0.036) |
| E. coli uniform | 6 | 0.032 (0.000) | 0.013 (0.006) | 0.129 (0.122) | 0.051 (0.044) | 0.013 (0.006) | 0.014 (0.006) | 0.013 (0.008) |
| E. coli uniform | 10 | 0.937 (0.000) | 0.424 (0.012) | 0.544 (0.128) | 0.481 (0.066) | 0.428 (0.013) | 0.423 (0.012) | 0.065 (0.024) |
| E. coli uniform | 15 | 39.430 (0.000) | 24.856 (0.026) | 24.779 (0.190) | 25.450 (0.120) | 25.065 (0.026) | 24.780 (0.025) | 3.861 (0.038) |
| stress clustered | 6 | 0.046 (0.000) | 0.023 (0.006) | 0.345 (0.328) | 0.080 (0.062) | 0.038 (0.021) | 0.224 (0.207) | 0.026 (0.012) |
| stress clustered | 10 | 1.095 (0.000) | 0.504 (0.013) | 1.199 (0.690) | 0.588 (0.091) | 0.529 (0.038) | 1.040 (0.546) | 0.097 (0.024) |
| stress clustered | 15 | 38.328 (0.000) | 25.669 (0.030) | 26.831 (0.833) | 26.673 (0.161) | 26.992 (0.028) | 26.290 (0.028) | 4.202 (0.039) |
