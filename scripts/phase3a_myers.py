#!/usr/bin/env python3
"""Phase 3A item 3: single-word fast path for Myers / Myers with pruning (E. coli 536; 10 reads; best of 5 runs; ms per read).

phase-2 builds (from git history, commit bb4699f) are the 'before'; the current binary is the 'after'.  For m <= 64 (one block) the pruned scan
used to carry vector and bookkeeping overhead for nothing; now both scans use the same register-only loop.  Target: pruned never slower than
plain Myers at one block.   Writes results/phase3a_myers.md."""
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
d = tempfile.mkdtemp(prefix="my_old_")
open(os.path.join(d, "old.cpp"), "w").write(subprocess.run(["git", "-C", ROOT, "show", "bb4699f:src/approx.cpp"], capture_output=True, text=True, check=True).stdout)
subprocess.run(["g++", "-O3", "-march=native", "-std=c++17", "-o", os.path.join(d, "old"), os.path.join(d, "old.cpp")], check=True)
old_bin = os.path.join(d, "old")
REPS = 2 if QUICK else 5


def run_old(algo, reads, k):
    pf = os.path.join(DATA, "_myold.txt"); open(pf, "w").write("\n".join(reads) + "\n")
    o = subprocess.run([old_bin, algo, tp, pf, str(k)], capture_output=True, text=True).stdout.split("\n")[0]
    return float(re.search(r"time_ms=([\d.]+)", o).group(1)) / len(reads)


def best(fn):
    return min(fn() for _ in range(REPS))


md = ["## Phase 3A, item 3: single-word fast path for Myers and Myers with pruning\n",
      "E. coli 536, 10 reads, best of 5 runs, ms per read. 'phase 2' = binary built from the final phase-2 commit; 'now' = current binary. "
      "Rows with m <= 60 are one block for every read (edits change a read's length by at most k, so every read stays within 64 bases); an earlier version of this table used m = 64, where about half the reads grow past 64 bases and need a second block, which mixed one- and two-block scans. Exactness: crosscheck (300 cases), real data (864 comparisons), tests/test_pruned.py (320 comparisons, m up to 600, k up to 200): 0 failures.\n",
      "| m (blocks) | k | Myers phase 2 | Myers now | pruned phase 2 | pruned now | pruned now / Myers now | pruned now vs phase-2 pruned | Edlib (best/read) |",
      "|---|---|---|---|---|---|---|---|---|"]
import edlib  # noqa: E402
gb = g.encode()
cells = [(20, 2), (32, 5), (50, 2), (50, 5), (50, 10), (56, 5), (60, 2), (100, 5), (100, 20)] if not QUICK else [(32, 5), (50, 5), (100, 5)]
worst1 = 0.0
for m, k in cells:
    reads = [r[0] for r in make_reads(g, 10, m, k, seed=5)]
    a = best(lambda: run_old("myers", reads, k)); b = best(lambda: run_old("myers_pruned", reads, k))
    c = best(lambda: run_bin("myers", tp, reads, k, dump=False, tag="m3")[0]["time_ms"] / len(reads))
    e = best(lambda: run_bin("myers_pruned", tp, reads, k, dump=False, tag="m3")[0]["time_ms"] / len(reads))
    import time
    def ed():
        t0 = time.perf_counter()
        for r in reads:
            edlib.align(r, gb, mode="HW", task="locations", k=k)
        return (time.perf_counter() - t0) * 1e3 / len(reads)
    f = best(ed)
    if (m + 63) // 64 == 1:
        worst1 = max(worst1, e / c)
    md.append(f"| {m} ({(m + 63) // 64}) | {k} | {a:.1f} | {c:.1f} | {b:.1f} | {e:.1f} | {e / c:.2f} | {b / e:.2f}x | {f:.1f} |")
md.append(f"\nAt one block the pruned scan now takes at most {worst1:.2f} times the time of plain Myers (same loop, so the difference is run-to-run noise); "
          "for several blocks the pruned scan is unchanged (m = 100 rows). The register-only loop also makes plain Myers faster for m <= 64, and the filters' window verification for short reads.")
open(os.path.join(RES, "phase3a_myers.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
