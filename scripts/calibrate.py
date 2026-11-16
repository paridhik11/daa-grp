#!/usr/bin/env python3
"""Fit the cost-model constants (ns) on the SYNTHETIC genome (data/genome_synth.txt).

The predicted-vs-measured evaluation uses E. coli 536, a different dataset (different text, different repeat
structure), so the evaluation is out-of-sample in text and read set (constants are still machine-specific).

  gen_ns   ~ c_hit * hits                      (hit enumeration + window sort/merge)
  ver_ns   ~ c_win * windows + c_bp * B * verified_bp
  scan     ~ a_scan * n * B                    (full Myers scan; B = ceil(m/64))
  sel_ns   ~ c_sel * m * 11                    (frequency-aware selection: counts table + DP)
Calibration k covers the filter/scan crossover (k up to 25 for m = 100): a first NNLS fit on k <= 15 only
underestimated verification cost near the crossover and the fallback misfired (see REPORT.md, Section 4).
Writes results/costmodel.json.  The constants are fitted in-sample on this machine (documented limitation).
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import nnls

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin  # noqa: E402
from simulate import make_reads  # noqa: E402

TEXT = os.path.join(DATA, "genome_synth.txt")
g = open(TEXT).read().strip()
n = len(g)
NP = 12
rows_gen, rows_ver, scan, sel = [], [], [], []
for m in (100, 200):
    B = (m + 63) // 64
    for k in (2, 4, 6, 8, 10, 12, 15, 18, 21, 25):
        reads = [r[0] for r in make_reads(g, NP, m, k, seed=777 + k + m)]
        for algo in ("filter", "filter_fa"):
            info, _ = run_bin(algo, TEXT, reads, k, nofallback=True, dump=False, tag="cal")
            rows_gen.append((info["hits"], info["gen_ms"] * 1e6))
            rows_ver.append((info["windows"], B * info["verified_bp"], info["ver_ms"] * 1e6))
            if algo == "filter_fa":
                sel.append((m * 11 * NP, info["sel_ms"] * 1e6))
    info, _ = run_bin("myers", TEXT, reads, 5, dump=False, tag="cal")
    scan.append((n * B * NP, info["time_ms"] * 1e6))
gh = np.array(rows_gen)
# all constants by non-negative least squares (c >= 0)
c_hit = float(nnls(gh[:, :1], gh[:, 1])[0][0])
A = np.array([[w, b] for w, b, _ in rows_ver]); y = np.array([t for *_, t in rows_ver])
sol, resid = nnls(A, y)
c_win, c_bp = float(sol[0]), float(sol[1])
# unconstrained least squares for comparison (documents whether the sign constraint is active)
sol_ls, *_ = np.linalg.lstsq(A, y, rcond=None)
print("unconstrained least squares c_win, c_bp =", sol_ls.tolist(), "; NNLS =", sol.tolist())
a_scan = float(sum(t for _, t in scan) / sum(x for x, _ in scan))
c_sel = float(sum(t for _, t in sel) / sum(x for x, _ in sel))
cm = dict(fit="NNLS (scipy.optimize.nnls)", c_win_unconstrained=float(sol_ls[0]), c_bp_unconstrained=float(sol_ls[1]), c_hit=c_hit, c_win=c_win, c_bp=c_bp, a_scan=a_scan, c_sel=c_sel, text_len=n,
          calibrated_on="synthetic genome_synth.txt", evaluated_on="E. coli 536")
json.dump(cm, open(os.path.join(RES, "costmodel.json"), "w"), indent=1)
print(json.dumps(cm, indent=1))
