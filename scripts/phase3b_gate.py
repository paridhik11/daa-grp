#!/usr/bin/env python3
"""Phase 3B item 5: gating the window-aware objective (filter_fa_w) so it only runs on repeat-heavy reads.

Gate: a read is repeat-heavy if the maximum hit count of its equal-split seeds is at least GATE_MIN and at least GATE_RATIO times the count expected
for a random seed (n / 4^s0).  Reports (1) how the gate splits the stress and ordinary reads, (2) timings before/after on the stress cases that failed
(pure filters, fallback off, mean ms per read, best of 3): filter, filter_fa (hits), filter_fa_w with the sort estimate and with the binned estimate,
filter_fa_wg (gated; binned), filter_auto.  Tuned on E. coli 536 and the synthetic stress text only.  Writes results/phase3b_gate.md."""
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import BIN, DATA, RES, cm_env, run_bin  # noqa: E402
from simulate import make_reads, make_reads_lowcomplex  # noqa: E402

QUICK = "--quick" in sys.argv
NR = 30 if QUICK else 100
ecoli = open(os.path.join(DATA, "ecoli536.txt")).read().strip(); ep = os.path.join(DATA, "ecoli536.txt")
st = open(os.path.join(DATA, "stress_text.txt")).read().strip(); sp = os.path.join(DATA, "stress_text.txt")
regions = [tuple(x) for x in json.load(open(os.path.join(DATA, "stress_regions.json")))]
md = ["## Phase 3B, item 5: gated window objective (filter_fa_wg)\n"]


def gate_stats(path, reads, k, ratio, gmin):
    pf = os.path.join(DATA, "_gate.txt"); open(pf, "w").write("\n".join(reads) + "\n")
    env = dict(os.environ, APPROX_DEBUG="2", APPROX_GATE=str(ratio), APPROX_GATEMIN=str(gmin), APPROX_CM=cm_env() or "")
    err = subprocess.run([BIN, "filter_fa_wg", path, pf, str(k)], capture_output=True, text=True, env=env).stderr.splitlines()
    g = [l for l in err if l.startswith("GATE")]
    vals = [(float(l.split("maxc=")[1].split()[0]), float(l.split("lam=")[1].split()[0]), int(l.split("heavy=")[1])) for l in g]
    return vals


def mean_ms(algo, path, reads, k, env=None, reps=3):
    for kk, v in (env or {}).items():
        os.environ[kk] = str(v)
    try:
        best = 1e18; extra = None
        for _ in range(reps):
            info, _s, tq = run_bin(algo, path, reads, k, nofallback=True, dump=False, tag="gt", times=True)
            m = float(np.mean([q["ns"] for q in tq])) / 1e6
            if m < best:
                best = m; extra = info
        return best, extra
    finally:
        for kk in (env or {}):
            os.environ.pop(kk, None)


scen = [("E. coli uniform", ep, lambda k: [r[0] for r in make_reads(ecoli, NR, 100, k, seed=11)]),
        ("stress uniform", sp, lambda k: [r[0] for r in make_reads(st, NR, 100, k, seed=31)]),
        ("stress clustered", sp, lambda k: [r[0] for r in make_reads(st, NR, 100, k, seed=32, mode="cluster")]),
        ("stress low-complexity edge", sp, lambda k: [r[0] for r in make_reads_lowcomplex(st, regions, NR, 100, k, "edge", seed=34)]),
        ("stress low-complexity inside", sp, lambda k: [r[0] for r in make_reads_lowcomplex(st, regions, NR, 100, k, "inside", seed=35)])]
md.append("### Gate: fraction of reads flagged repeat-heavy (GATE_RATIO = 20, GATE_MIN = 100) and the ratio max seed count / random expectation (median, 90th percentile)\n")
md.append("| scenario | k | reads flagged | ratio median | ratio p90 |")
md.append("|---|---|---|---|---|")
for name, path, mk in scen:
    for k in (6, 10, 15):
        v = gate_stats(path, mk(k), k, 20, 100)
        ratios = np.array([mc / max(1.0, lam) for mc, lam, _ in v])
        md.append(f"| {name} | {k} | {sum(h for *_, h in v)} of {len(v)} | {np.median(ratios):.1f} | {np.percentile(ratios, 90):.1f} |")
md.append("\n### Before and after (pure filters, mean ms per read, best of 3; selection ms in brackets)\n")
VARS = [("filter (equal split)", "filter", None), ("filter_fa (hits)", "filter_fa", None),
        ("filter_fa_w, sort estimate (phase 2)", "filter_fa_w", {"APPROX_WBIN": 0}), ("filter_fa_w, binned estimate", "filter_fa_w", {"APPROX_WBIN": 1}),
        ("filter_fa_wg (gated, binned)", "filter_fa_wg", {"APPROX_WBIN": 1}), ("filter_fa_wg (gated, sort estimate)", "filter_fa_wg", {"APPROX_WBIN": 0}),
        ("filter_auto", "filter_auto", None)]
md.append("| scenario | k | " + " | ".join(v[0] for v in VARS) + " |")
md.append("|---|---|" + "---|" * len(VARS))
for name, path, mk in scen[3:] + scen[:1] + scen[2:3]:
    for k in ((10,) if QUICK else (6, 10, 15)):
        reads = mk(k); cells = []
        for vn, algo, env in VARS:
            m, info = mean_ms(algo, path, reads, k, env)
            cells.append(f"{m:.3f} ({info['sel_ms'] / len(reads):.3f})")
        md.append(f"| {name} | {k} | " + " | ".join(cells) + " |")
open(os.path.join(RES, "phase3b_gate.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
