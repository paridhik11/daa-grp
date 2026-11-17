## Phase 3A, item 2: faster Landau-Vishkin

E. coli 536, m = 100, 5 reads, best of 3 runs, ms per read. Every variant returns the same matches as Sellers (tests/test_lv.py: 600 comparisons on tandem-repeat, low-complexity and N texts; crosscheck: 300 random cases; real data: 864 comparisons; the scalar path is also tested with APPROX_LVSIMD=0).

| k | phase-2 LV | lv_direct (control) | lv_stream | lv scalar blocked | lv AVX-512 (default) | speedup of default vs phase-2 | Myers |
|---|---|---|---|---|---|---|---|
| 1 | 92.8 | 90.3 | 53.5 | 40.9 | 15.5 | 5.97x | 61.2 |
| 5 | 212.3 | 214.4 | 148.4 | 121.6 | 40.5 | 5.24x | 61.5 |
| 10 | 426.8 | 425.8 | 270.6 | 237.5 | 75.4 | 5.66x | 67.2 |
| 30 | 1189.8 | 1175.7 | 699.3 | 706.4 | 207.9 | 5.72x | 69.4 |

Variants and what they showed: the packed XOR + ctz extension alone (lv_stream) gave about 1.5x; temporal blocking (a chunk of diagonals carried through all layers in L1) gave a small and inconsistent gain (the scalar cell is instruction-bound, about 30 instructions or 4.4 ns, not memory-bound); the AVX-512 cell loop (8 cells per iteration, two gathers, lzcnt) gave the large step. A CPU without AVX-512 uses the scalar blocked loop (about 1.7x to 2.3x faster than phase 2). Texts or patterns containing N use the scalar loop.
