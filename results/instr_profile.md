## Instructions per cell and branch behaviour: banded vs Myers

Hardware performance counters: perf installed = False, CPU PMU device present = False. In this VM neither is true, so hardware counters could not be used. The numbers below come from `valgrind --tool=cachegrind --branch-sim=yes`: instruction counts are exact executed-instruction counts; branch mispredictions come from cachegrind's simulated predictor and are indicative only. They do not show cache misses or real timing.

Text: first 150,000 bp of E. coli 536; m = 100; per pattern (mean of 5). Myers cells = m * n DP cells represented (64 per word operation); banded cells = cells actually computed.

| k | algorithm | instructions per text char | instructions per cell | cells per text char | branches per text char | simulated mispredicts per 1000 chars | simulated mispredict rate |
|---|---|---|---|---|---|---|---|
| 5 | banded | 229.8 | 21.43 | 10.7 | 16.29 | 2024.7 | 12.43% |
| 5 | myers | 98.7 | 0.99 | 100.0 | 9.01 | 621.8 | 6.90% |
| 20 | banded | 754.3 | 18.92 | 39.9 | 45.42 | 1952.4 | 4.30% |
| 20 | myers | 98.8 | 0.99 | 100.0 | 9.01 | 639.4 | 7.10% |

At k = 5: banded executes 2.3x the instructions of Myers per text character while computing 9x fewer cells; instructions per cell are 21.4 (banded) against 0.99 (Myers).

At k = 20: banded executes 7.6x the instructions of Myers per text character while computing 3x fewer cells; instructions per cell are 18.9 (banded) against 0.99 (Myers).
