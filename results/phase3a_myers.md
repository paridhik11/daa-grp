## Phase 3A, item 3: single-word fast path for Myers and Myers with pruning

E. coli 536, 10 reads, best of 5 runs, ms per read. 'phase 2' = binary built from the final phase-2 commit; 'now' = current binary. Rows with m <= 60 are one block for every read (edits change a read's length by at most k, so every read stays within 64 bases); an earlier version of this table used m = 64, where about half the reads grow past 64 bases and need a second block, which mixed one- and two-block scans. Exactness: crosscheck (300 cases), real data (864 comparisons), tests/test_pruned.py (320 comparisons, m up to 600, k up to 200): 0 failures.

| m (blocks) | k | Myers phase 2 | Myers now | pruned phase 2 | pruned now | pruned now / Myers now | pruned now vs phase-2 pruned | Edlib (best/read) |
|---|---|---|---|---|---|---|---|---|
| 20 (1) | 2 | 39.6 | 19.3 | 46.3 | 18.6 | 0.96 | 2.49x | 43.8 |
| 32 (1) | 5 | 38.9 | 19.0 | 48.6 | 19.4 | 1.02 | 2.50x | 34.4 |
| 50 (1) | 2 | 39.7 | 19.4 | 48.5 | 19.6 | 1.01 | 2.47x | 34.6 |
| 50 (1) | 5 | 37.6 | 18.7 | 50.2 | 18.9 | 1.02 | 2.65x | 35.4 |
| 50 (1) | 10 | 40.8 | 19.6 | 48.8 | 19.0 | 0.97 | 2.58x | 33.2 |
| 56 (1) | 5 | 39.0 | 19.2 | 45.8 | 18.7 | 0.98 | 2.44x | 34.7 |
| 60 (1) | 2 | 37.7 | 18.7 | 49.2 | 18.8 | 1.00 | 2.62x | 35.9 |
| 100 (2) | 5 | 66.1 | 68.8 | 46.9 | 46.4 | 0.67 | 1.01x | 34.8 |
| 100 (2) | 20 | 67.9 | 75.0 | 46.9 | 48.6 | 0.65 | 0.96x | 36.2 |

At one block the pruned scan now takes at most 1.02 times the time of plain Myers (same loop, so the difference is run-to-run noise); for several blocks the pruned scan is unchanged (m = 100 rows). The register-only loop also makes plain Myers faster for m <= 64, and the filters' window verification for short reads.
