#!/usr/bin/env python3
"""Phase 3B: does the cost model support adding the s-hit filters to filter_auto?  Fallback active everywhere (real use); mean ms per read, best of 3.
filter_auto(no s) = APPROX_AUTOS=0 (the phase-3A auto), filter_auto = with filter_fs2 / filter_fs3 as options chosen by the predicted cost
(Poisson-binomial window prediction x 2.0, tuned on E. coli 536 only).  `auto/best` = filter_auto / min(filter, filter_fa, filter_fs2, filter_fs3, Myers).
Writes results/phase3b_auto.md."""
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin  # noqa: E402
from simulate import make_reads  # noqa: E402

QUICK = "--quick" in sys.argv
NR = 30 if QUICK else 100
REPS = 1 if QUICK else 3
md = ["## Phase 3B: s-hit filters inside filter_auto\n",
      "Mean ms per read (fallback active), best of %d runs, %d reads. Choices: queries that filter_auto sent to the 2-hit / 3-hit filter (of all queries).\n" % (REPS, NR),
      "| text | setting | Myers | filter | filter_fa | filter_fs2 | filter_fs3 | filter_auto (no s) | filter_auto | auto / best fixed | choices s2 / s3 |", "|---|---|---|---|---|---|---|---|---|---|---|"]
worst = []


def mean(algo, path, reads, k, env=None):
    for kk, v in (env or {}).items():
        os.environ[kk] = str(v)
    best, info = 1e18, None
    for _ in range(REPS):
        i, _s, tq = run_bin(algo, path, reads, k, dump=False, tag="pa", times=True)
        m = float(np.mean([q["ns"] for q in tq])) / 1e6
        if m < best:
            best, info = m, i
    for kk in (env or {}):
        os.environ.pop(kk, None)
    return best, info


def row(label, path, reads, k, setting):
    t = {}
    for name, algo, env in (("my", "myers", None), ("eq", "filter", None), ("fa", "filter_fa", None), ("s2", "filter_fs2", None), ("s3", "filter_fs3", None),
                            ("a0", "filter_auto", {"APPROX_AUTOS": 0}), ("a1", "filter_auto", None)):
        t[name] = mean(algo, path, reads, k, env)
    best = min(t[x][0] for x in ("my", "eq", "fa", "s2", "s3"))
    ratio = t["a1"][0] / best
    worst.append((ratio, label, setting))
    info = t["a1"][1]
    md.append(f"| {label} | {setting} | {t['my'][0]:.3f} | {t['eq'][0]:.3f} | {t['fa'][0]:.3f} | {t['s2'][0]:.3f} | {t['s3'][0]:.3f} | {t['a0'][0]:.3f} | {t['a1'][0]:.3f} | {ratio:.2f} | {info.get('auto_s2', 0):.0f} / {info.get('auto_s3', 0):.0f} |")


ecoli = open(os.path.join(DATA, "ecoli536.txt")).read().strip(); ep = os.path.join(DATA, "ecoli536.txt")
for k in ([8, 12, 15] if QUICK else [6, 8, 10, 12, 14, 15, 17, 20]):
    row("E. coli 536", ep, [r[0] for r in make_reads(ecoli, NR, 100, k, seed=11)], k, f"m=100, k={k}")
for m, k in ([(200, 30)] if QUICK else [(200, 20), (200, 30), (400, 60)]):
    row("E. coli 536", ep, [r[0] for r in make_reads(ecoli, NR, m, k, seed=12)], k, f"m={m}, k={k}")
sp = os.path.join(DATA, "genome_synth.txt")
synth = open(sp).read().strip()
for k in ([12] if QUICK else [10, 12, 15]):
    row("synthetic", sp, [r[0] for r in make_reads(synth, NR, 100, k, seed=21)], k, f"m=100, k={k}")
human = open(os.path.join(DATA, "mt_human.txt")).read().strip(); orang = open(os.path.join(DATA, "mt_orang.txt")).read().strip(); hp = os.path.join(DATA, "mt_human.txt")
for m, ks in ((100, [15, 20]), (200, [30])):
    rng = random.Random(100 + m)
    segs = [orang[s:s + m] for s in (rng.randrange(0, len(orang) - m) for _ in range(NR))]
    for k in ks:
        row("mt_orang vs mt_human (real)", hp, segs, k, f"m={m}, k={k}")
worst.sort(reverse=True)
md.append(f"\nWorst auto/best-fixed ratio: {worst[0][0]:.2f} ({worst[0][1]}, {worst[0][2]}); {sum(1 for w in worst if w[0] > 1.05)} of {len(worst)} settings exceed 1.05.")
open(os.path.join(RES, "phase3b_auto.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
