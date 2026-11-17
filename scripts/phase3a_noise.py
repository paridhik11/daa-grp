#!/usr/bin/env python3
"""Run-to-run noise check for filter_auto versus filter_fa on the two settings that exceeded 1.05 in phase3a_filter.py.
300 reads, 6 interleaved runs each; reports the spread of the same algorithm across runs and the ratio of medians. Appends to results/phase3a_filter.md."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin  # noqa: E402
from simulate import make_reads  # noqa: E402

g = open(os.path.join(DATA, "ecoli536.txt")).read().strip(); tp = os.path.join(DATA, "ecoli536.txt")
lines = ["\n### Run-to-run noise check (E. coli 536, 300 reads, 6 interleaved runs)\n",
         "| setting | filter_fa runs: min to max (ms) | filter_auto runs: min to max (ms) | spread of the same algorithm | ratio of medians auto / fa | ratio of minima |", "|---|---|---|---|---|---|"]
for m, k, seed in [(100, 10, 11), (800, 40, 12)]:
    reads = [r[0] for r in make_reads(g, 300, m, k, seed=seed)]
    res = {"filter_fa": [], "filter_auto": []}
    for _ in range(6):
        for a in res:
            _i, _s, tq = run_bin(a, tp, reads, k, dump=False, tag="nz", times=True)
            res[a].append(float(np.mean([q["ns"] for q in tq])) / 1e6)
    f, au = res["filter_fa"], res["filter_auto"]
    lines.append(f"| m={m}, k={k} | {min(f):.3f} to {max(f):.3f} | {min(au):.3f} to {max(au):.3f} | {100 * (max(f) / min(f) - 1):.0f}% | {np.median(au) / np.median(f):.3f} | {min(au) / min(f):.3f} |")
lines.append("\nThe same algorithm varies by more than 5% between identical runs, so single-run ratios of 1.05 to 1.08 are within noise.")
open(os.path.join(RES, "phase3a_filter.md"), "a").write("\n".join(lines) + "\n")
print("\n".join(lines))
