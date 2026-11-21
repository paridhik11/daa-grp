#!/usr/bin/env python3
"""Phase 3B item 4: generalised pigeonhole (s required seed hits in a diagonal band) on E. coli 536 reads at 10 to 25% error and on real
homologous mt_orang segments against mt_human.  Pure filters (fallback off) so candidate counts are visible; exactness asserted against Myers.
Tuned/evaluated on E. coli 536 and mt data only (never B. fragilis).   Writes results/phase3b_filter_s.md and .csv.   Usage: [--quick]"""
import csv
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin  # noqa: E402
from simulate import make_reads  # noqa: E402

QUICK = "--quick" in sys.argv
ALGOS = [("Myers scan", "myers"), ("filter (equal, s=1)", "filter"), ("filter_fa (fa, s=1)", "filter_fa"), ("filter_s2 (equal, s=2)", "filter_s2"),
         ("filter_s3 (equal, s=3)", "filter_s3"), ("filter_fs2 (fa, s=2)", "filter_fs2"), ("filter_fs3 (fa, s=3)", "filter_fs3")]
rows, md = [], ["## Phase 3B, item 4: generalised pigeonhole (s required hits in a diagonal band of width k)\n",
                "k + s disjoint seeds; a candidate needs s distinct seeds hit within a band of k diagonals; s = 1 is the classic filter. Pure filters (cost-model fallback off), "
                "mean per read; `windows` = merged candidate windows verified, `bp` = bases verified, `gen` = time to build candidates (sort and band pass). "
                "Every filter's matches are asserted equal to Myers on every read.\n"]


def run_setting(label, tpath, reads, k):
    md.append(f"\n### {label}, k = {k}\n")
    md.append("| method | ms per read | seed hits | windows | bp verified | gen ms | selection ms |")
    md.append("|---|---|---|---|---|---|---|")
    truth = None
    for name, algo in ALGOS:
        info, sets = run_bin(algo, tpath, reads, k, nofallback=True, tag="p3b")
        if algo == "myers":
            truth = sets
        elif {i: v for i, v in sets.items() if v} != {i: v for i, v in truth.items() if v}:
            raise SystemExit(f"EXACT MISMATCH {label} k={k} {algo}")
        nr = len(reads)
        # time without --dump overhead of the sets: rerun with times
        best = min(run_bin(algo, tpath, reads, k, nofallback=True, dump=False, tag="p3b")[0]["time_ms"] for _ in range(2)) / nr
        r = dict(setting=label, k=k, method=name, ms=best, hits=info["hits"] / nr, windows=info["windows"] / nr, bp=info["verified_bp"] / nr,
                 gen_ms=info["gen_ms"] / nr, sel_ms=info["sel_ms"] / nr)
        rows.append(r)
        md.append(f"| {name} | {best:.3f} | {r['hits']:.0f} | {r['windows']:.0f} | {r['bp']:.0f} | {r['gen_ms']:.3f} | {r['sel_ms']:.3f} |")


ecoli = open(os.path.join(DATA, "ecoli536.txt")).read().strip(); ep = os.path.join(DATA, "ecoli536.txt")
NR = 20 if QUICK else 50
for k in ([12, 20] if QUICK else [10, 12, 15, 17, 20, 25]):
    run_setting("E. coli 536, m = 100 (%d%% error)" % k, ep, [r[0] for r in make_reads(ecoli, NR, 100, k, seed=11)], k)
human = open(os.path.join(DATA, "mt_human.txt")).read().strip(); orang = open(os.path.join(DATA, "mt_orang.txt")).read().strip()
hp = os.path.join(DATA, "mt_human.txt")
for m, ks in ((100, [15, 20, 25, 30]), (200, [30, 40])):
    rng = random.Random(100 + m)
    segs = [orang[s:s + m] for s in (rng.randrange(0, len(orang) - m) for _ in range(100 if not QUICK else 30))]
    for k in (ks[:2] if QUICK else ks):
        run_setting(f"mt_orang segments vs mt_human (real homologs), m = {m}", hp, segs, k)
open(os.path.join(RES, "phase3b_filter_s.md"), "w").write("\n".join(md) + "\n")
with open(os.path.join(RES, "phase3b_filter_s.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("\n".join(md))
