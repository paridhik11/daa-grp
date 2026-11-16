#!/usr/bin/env python3
"""Instructions per DP cell and branch mispredictions for banded vs Myers.

perf hardware counters are unavailable in this VM (no PMU device under /sys/bus/event_source/devices), so this uses
valgrind --tool=cachegrind --branch-sim=yes on a generic x86-64 build (valgrind cannot run -march=native AVX-512 code):
  * instruction counts (Ir) are exact counts of executed x86 instructions;
  * branch mispredictions come from cachegrind's SIMULATED predictor, not from hardware (indicative only).
Per-pattern cost = (run with 6 patterns - run with 1 pattern) / 5, which removes text loading and start-up.
Writes results/instr_profile.md (appended to summary.md by plot.py).
"""
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, ROOT  # noqa: E402
from simulate import make_reads  # noqa: E402

lines = ["## Instructions per cell and branch behaviour: banded vs Myers\n"]
pmu = os.path.isdir("/sys/bus/event_source/devices/cpu") or os.path.isdir("/sys/bus/event_source/devices/cpu_core")
lines.append(f"Hardware performance counters: perf installed = {bool(shutil.which('perf'))}, CPU PMU device present = {pmu}. "
             "In this VM neither is true, so hardware counters could not be used. The numbers below come from `valgrind --tool=cachegrind "
             "--branch-sim=yes`: instruction counts are exact executed-instruction counts; branch mispredictions come from cachegrind's "
             "simulated predictor and are indicative only. They do not show cache misses or real timing.\n")
if not shutil.which("valgrind"):
    lines.append("valgrind is not installed: instruction-level measurement skipped. The Myers versus banded explanation stays a hypothesis.\n")
    open(os.path.join(RES, "instr_profile.md"), "w").write("\n".join(lines) + "\n")
    sys.exit(0)

gen = os.path.join(ROOT, "approx_generic")
subprocess.run(["g++", "-O3", "-std=c++17", "-o", gen, os.path.join(ROOT, "src", "approx.cpp")], check=True)
g = open(os.path.join(DATA, "genome.txt")).read().strip()
N = 150_000
text = os.path.join(DATA, "_valgrind_text.txt")
open(text, "w").write(g[:N] + "\n")
m = 100


def cg(algo, pats, k):
    pf = os.path.join(DATA, "_valgrind_pats.txt")
    open(pf, "w").write("\n".join(pats) + "\n")
    r = subprocess.run(["valgrind", "--tool=cachegrind", "--cache-sim=no", "--branch-sim=yes", "--cachegrind-out-file=/dev/null",
                        gen, algo, text, pf, str(k)], capture_output=True, text=True)
    err = r.stderr
    num = lambda pat: int(re.search(pat, err).group(1).replace(",", ""))
    ir = num(r"I\s+refs:\s+([\d,]+)")
    br = num(r"Branches:\s+([\d,]+)")
    mis = num(r"Mispredicts:\s+([\d,]+)")
    cells = re.search(r"cells_banded=(\d+)", r.stdout)
    return ir, br, mis, int(cells.group(1)) if cells else 0


rows = []
for k in (5, 20):
    reads = [r[0] for r in make_reads(g[:N], 6, m, k, seed=99)]
    res = {}
    for algo in ("banded", "myers"):
        a6, a1 = cg(algo, reads, k), cg(algo, reads[:1], k)
        d = [(x - y) / 5.0 for x, y in zip(a6[:3], a1[:3])]
        cells = (a6[3] - a1[3]) / 5.0 if algo == "banded" else float(m * N)
        res[algo] = dict(ir=d[0], br=d[1], mis=d[2], cells=cells)
    rows.append((k, res))

lines.append(f"Text: first {N:,} bp of E. coli 536; m = {m}; per pattern (mean of 5). Myers cells = m * n DP cells represented "
             "(64 per word operation); banded cells = cells actually computed.\n")
lines.append("| k | algorithm | instructions per text char | instructions per cell | cells per text char | branches per text char | simulated mispredicts per 1000 chars | simulated mispredict rate |")
lines.append("|---|---|---|---|---|---|---|---|")
for k, res in rows:
    for algo in ("banded", "myers"):
        r = res[algo]
        lines.append(f"| {k} | {algo} | {r['ir'] / N:.1f} | {r['ir'] / r['cells']:.2f} | {r['cells'] / N:.1f} | {r['br'] / N:.2f} | "
                     f"{1000 * r['mis'] / N:.1f} | {100 * r['mis'] / max(1, r['br']):.2f}% |")
for k, res in rows:
    b, my = res["banded"], res["myers"]
    lines.append(f"\nAt k = {k}: banded executes {b['ir'] / my['ir']:.1f}x the instructions of Myers per text character "
                 f"while computing {my['cells'] / b['cells']:.0f}x fewer cells; instructions per cell are {b['ir'] / b['cells']:.1f} (banded) "
                 f"against {my['ir'] / my['cells']:.2f} (Myers).")
open(os.path.join(RES, "instr_profile.md"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
