#!/usr/bin/env python3
"""Final pass: the repeat gate of filter_fa_wg on B. fragilis (never used for tuning), thresholds frozen (ratio 20, min count 100).
Counts reads flagged repeat-heavy and the ratio max equal-split seed count / random expectation. Writes results/final_gate.md."""
import os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, BIN, cm_env  # noqa: E402
from simulate import make_reads  # noqa: E402
path = os.path.join(DATA, "bfragilis.txt"); t = open(path).read().strip()
md = ["## Final pass: gate on B. fragilis (BLIND, thresholds frozen: ratio 20, minimum count 100)\n", "| k | reads flagged | ratio median | ratio p90 | ratio max |", "|---|---|---|---|---|"]
for k in (6, 10, 15):
    reads = [r[0] for r in make_reads(t, 200, 100, k, seed=303 + k)]
    pf = os.path.join(DATA, "_fg.txt"); open(pf, "w").write("\n".join(reads) + "\n")
    env = dict(os.environ, APPROX_DEBUG="2", APPROX_CM=cm_env() or "")
    err = subprocess.run([BIN, "filter_fa_wg", path, pf, str(k)], capture_output=True, text=True, env=env).stderr.splitlines()
    v = [(float(l.split("maxc=")[1].split()[0]), float(l.split("lam=")[1].split()[0]), int(l.split("heavy=")[1])) for l in err if l.startswith("GATE")]
    r = np.array([a / max(1.0, b) for a, b, _ in v])
    md.append(f"| {k} | {sum(h for *_, h in v)} of {len(v)} | {np.median(r):.1f} | {np.percentile(r, 90):.1f} | {r.max():.1f} |")
open(os.path.join(RES, "final_gate.md"), "w").write("\n".join(md) + "\n"); print("\n".join(md))
