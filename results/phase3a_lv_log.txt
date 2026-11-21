## Phase 3A, item 2: faster Landau-Vishkin

E. coli 536, m = 100, 5 reads, best of 3 runs, ms per read. Every variant returns the same matches as Sellers (tests/test_lv.py: 600 comparisons on tandem-repeat, low-complexity and N texts; crosscheck: 300 random cases; real data: 864 comparisons; the scalar path is also tested with APPROX_LVSIMD=0).

| k | phase-2 LV | lv_direct (control) | lv_stream | lv scalar blocked | lv AVX-512 (default) | speedup of default vs phase-2 | Myers |
|---|---|---|---|---|---|---|---|
| 1 | 94.1 | 94.6 | 66.7 | 46.5 | 16.5 | 5.70x | 70.2 |
| 5 | 239.6 | 240.0 | 165.8 | 135.2 | 43.6 | 5.49x | 70.0 |
| 10 | 444.0 | 443.9 | 280.6 | 244.6 | 79.9 | 5.56x | 73.1 |
| 30 | 1232.2 | 1245.9 | 758.2 | 695.0 | 221.5 | 5.56x | 72.5 |

Variants and what they showed: the packed XOR + ctz extension alone (lv_stream) gave about 1.5x; temporal blocking (a chunk of diagonals carried through all layers in L1) gave a small and inconsistent gain (the scalar cell is instruction-bound, about 30 instructions or 4.4 ns, not memory-bound); the AVX-512 cell loop (8 cells per iteration, two gathers, lzcnt) gave the large step. A CPU without AVX-512 uses the scalar blocked loop (about 1.7x to 2.3x faster than phase 2). Texts or patterns containing N use the scalar loop.
