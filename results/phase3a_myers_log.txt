## Phase 3A, item 3: single-word fast path for Myers and Myers with pruning

E. coli 536, 10 reads, best of 5 runs, ms per read. 'phase 2' = binary built from the final phase-2 commit; 'now' = current binary. Rows with m <= 60 are one block for every read (edits change a read's length by at most k, so every read stays within 64 bases); an earlier version of this table used m = 64, where about half the reads grow past 64 bases and need a second block, which mixed one- and two-block scans. Exactness: crosscheck (300 cases), real data (864 comparisons), tests/test_pruned.py (320 comparisons, m up to 600, k up to 200): 0 failures.

| m (blocks) | k | Myers phase 2 | Myers now | pruned phase 2 | pruned now | pruned now / Myers now | pruned now vs phase-2 pruned | Edlib (best/read) |
|---|---|---|---|---|---|---|---|---|
| 20 (1) | 2 | 41.2 | 20.6 | 50.2 | 20.6 | 1.00 | 2.43x | 36.5 |
| 32 (1) | 5 | 41.1 | 20.6 | 50.9 | 20.6 | 1.00 | 2.47x | 38.2 |
| 50 (1) | 2 | 41.1 | 20.7 | 50.1 | 20.6 | 0.99 | 2.43x | 37.7 |
| 50 (1) | 5 | 41.7 | 20.6 | 51.7 | 20.6 | 1.00 | 2.51x | 36.9 |
| 50 (1) | 10 | 41.3 | 20.6 | 50.4 | 20.6 | 1.00 | 2.44x | 37.9 |
| 56 (1) | 5 | 41.4 | 20.8 | 50.8 | 20.8 | 1.00 | 2.44x | 36.2 |
| 60 (1) | 2 | 41.5 | 20.8 | 51.8 | 20.8 | 1.00 | 2.48x | 36.7 |
| 100 (2) | 5 | 69.3 | 72.5 | 49.8 | 52.2 | 0.72 | 0.95x | 38.6 |
| 100 (2) | 20 | 68.6 | 71.6 | 50.2 | 53.0 | 0.74 | 0.95x | 38.6 |

At one block the pruned scan now takes at most 1.00 times the time of plain Myers (same loop, so the difference is run-to-run noise); for several blocks the pruned scan is unchanged (m = 100 rows). The register-only loop also makes plain Myers faster for m <= 64, and the filters' window verification for short reads.
