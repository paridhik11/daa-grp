#!/usr/bin/env python3
"""Phase 3A item 2: Landau-Vishkin variants (E. coli 536, m = 100, 5 reads, best of 3 runs, ms per read).

  phase-2 LV : the implementation from the final phase-2 commit, built from git history (direct character extension, streamed layers)
  lv_direct  : the same code in the current binary (control)
  lv_stream  : 2-bit packed text/pattern, XOR + ctz extension, three-way max in a vectorisable pass, layers streamed over the text
  lv scalar  : + temporal blocking (2048-diagonal chunks through all layers in L1 buffers), scalar cell loop (APPROX_LVSIMD=0)
  lv AVX-512 : + 8 cells per iteration with gathers and lzcnt (the default `lv`)
Exactness of every variant is checked against Sellers by tests/test_lv.py, tests/crosscheck.py and tests/crosscheck_real.py.
Writes results/phase3a_lv.md."""
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, ROOT, run_bin  # noqa: E402
from simulate import make_reads  # noqa: E402

QUICK = "--quick" in sys.argv
g = open(os.path.join(DATA, "ecoli536.txt")).read().strip(); tp = os.path.join(DATA, "ecoli536.txt")
old_bin = None
try:
    d = tempfile.mkdtemp(prefix="lv_old_")
    src = subprocess.run(["git", "-C", ROOT, "show", "bb4699f:src/approx.cpp"], capture_output=True, text=True, check=True).stdout
    open(os.path.join(d, "old.cpp"), "w").write(src)
    subprocess.run(["g++", "-O3", "-march=native", "-std=c++17", "-o", os.path.join(d, "old"), os.path.join(d, "old.cpp")], check=True)
    old_bin = os.path.join(d, "old")
except Exception as ex:  # noqa: BLE001
    print("could not build the phase-2 binary:", ex)


def run_old(reads, k):
    pf = os.path.join(DATA, "_lvold.txt"); open(pf, "w").write("\n".join(reads) + "\n")
    o = subprocess.run([old_bin, "lv", tp, pf, str(k)], capture_output=True, text=True).stdout.split("\n")[0]
    return float(re.search(r"time_ms=([\d.]+)", o).group(1)) / len(reads)


def best(fn, reps=3):
    return min(fn() for _ in range(reps))


md = ["## Phase 3A, item 2: faster Landau-Vishkin\n",
      "E. coli 536, m = 100, 5 reads, best of 3 runs, ms per read. Every variant returns the same matches as Sellers (tests/test_lv.py: 600 comparisons on tandem-repeat, low-complexity and N texts; crosscheck: 300 random cases; real data: 864 comparisons; the scalar path is also tested with APPROX_LVSIMD=0).\n",
      "| k | phase-2 LV | lv_direct (control) | lv_stream | lv scalar blocked | lv AVX-512 (default) | speedup of default vs phase-2 | Myers |", "|---|---|---|---|---|---|---|---|"]
for k in ((1, 5) if QUICK else (1, 5, 10, 30)):
    reads = [r[0] for r in make_reads(g, 5, 100, k, seed=3)]
    old = best(lambda: run_old(reads, k)) if old_bin else float("nan")
    t = {}
    for name, algo, env in (("direct", "lv_direct", None), ("stream", "lv_stream", None), ("scalar", "lv", "0"), ("simd", "lv", None), ("myers", "myers", None)):
        if env is not None:
            os.environ["APPROX_LVSIMD"] = env
        t[name] = best(lambda: run_bin(algo, tp, reads, k, dump=False, tag="lv3")[0]["time_ms"] / len(reads))
        os.environ.pop("APPROX_LVSIMD", None)
    md.append(f"| {k} | {old:.1f} | {t['direct']:.1f} | {t['stream']:.1f} | {t['scalar']:.1f} | {t['simd']:.1f} | {old / t['simd']:.2f}x | {t['myers']:.1f} |")
md.append("\nVariants and what they showed: the packed XOR + ctz extension alone (lv_stream) gave about 1.5x; temporal blocking (a chunk of diagonals carried through all layers in L1) "
          "gave a small and inconsistent gain (the scalar cell is instruction-bound, about 30 instructions or 4.4 ns, not memory-bound); the AVX-512 cell loop (8 cells per iteration, "
          "two gathers, lzcnt) gave the large step. A CPU without AVX-512 uses the scalar blocked loop (about 1.7x to 2.3x faster than phase 2). Texts or patterns containing N use the scalar loop.")
open(os.path.join(RES, "phase3a_lv.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
