#!/usr/bin/env python3
"""Does changing the seed-selection objective from `hits` to `hits + lambda * windows` help?  (item 5)

filter     : equal split
filter_fa  : DP minimises total seed hits (before)
filter_fa_w: DP minimises hits + lambda * windows, lambda = (c_win + c_bp B (m + 2k)) / c_hit from the fitted cost model (after)
All pure filters (fallback off), exact (asserted equal to Myers).  Stress scenarios on E. coli 536 + injected low-complexity arrays
(same reads as the main stress benchmark), plus an unmodified E. coli 536 k sweep with 100 reads.
Writes results/objective_test.md and results/objective_test.csv.
"""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin, load_cm  # noqa: E402
from simulate import inject_low_complexity, make_reads, make_reads_lowcomplex  # noqa: E402

ALGOS = [("filter", "filter", "equal split", None), ("filter_fa", "filter_fa", "fa: min hits", None),
         ("filter_fa_w", "filter_fa_w", "fa_w: min hits + lambda*windows (hit-list cap 20000)", None),
         ("filter_fa_w1k", "filter_fa_w", "fa_w, hit-list cap 1000 (tried after seeing the first result)", "1000")]
CM = load_cm()
main = open(os.path.join(DATA, "ecoli536.txt")).read().strip()
sp = os.path.join(DATA, "stress_text.txt")
rp = os.path.join(DATA, "stress_regions.json")
if not (os.path.exists(sp) and os.path.exists(rp)):
    st, regions = inject_low_complexity(main, 0.01, seed=5)
    open(sp, "w").write(st + "\n"); json.dump(regions, open(rp, "w"))
st = open(sp).read().strip()
regions = [tuple(x) for x in json.load(open(rp))]
rows, md = [], ["## Seed-selection objective: hits versus hits + lambda * windows\n"]
md.append("Pure filters (fallback off), mean over reads; windows = merged windows verified; lambda = (c_win + c_bp B (m + 2k)) / c_hit "
          f"with the fitted constants (c_hit = {CM['c_hit']:.1f}, c_win = {CM['c_win']:.0f}, c_bp = {CM['c_bp']:.2f} ns), e.g. lambda = "
          f"{(CM['c_win'] + CM['c_bp'] * 2 * (100 + 20)) / CM['c_hit']:.1f} at m = 100, k = 10. All methods return the same matches as Myers (asserted).\n")


def run_case(label, tpath, reads, k):
    _, truth = run_bin("myers", tpath, reads, k, tag="ot")
    out = {}
    for key, algo, desc, cap in ALGOS:
        if cap:
            os.environ["APPROX_WCAP"] = cap
        info, sets = run_bin(algo, tpath, reads, k, nofallback=True, tag="ot")
        os.environ.pop("APPROX_WCAP", None)
        if {i: s for i, s in sets.items() if s} != {i: s for i, s in truth.items() if s}:
            raise SystemExit(f"EXACT MISMATCH {label} k={k} {algo}")
        n = len(reads)
        out[key] = dict(hits=info["hits"] / n, windows=info["windows"] / n, bp=info["verified_bp"] / n, ms=info["time_ms"] / n, sel=info["sel_ms"] / n)
        rows.append(dict(case=label, k=k, algo=key, **out[key]))
    return out


md.append("### Stress scenarios (E. coli 536 + about 1% microsatellite/homopolymer arrays, m = 100, 20 reads)\n")
md.append("| scenario | k | algorithm | hits | windows | bp verified | query ms | seed selection ms |")
md.append("|---|---|---|---|---|---|---|---|")
scen = {"uniform": lambda k: [r[0] for r in make_reads(st, 20, 100, k, seed=31)],
        "clustered": lambda k: [r[0] for r in make_reads(st, 20, 100, k, seed=32, mode="cluster")],
        "indel_heavy": lambda k: [r[0] for r in make_reads(st, 20, 100, k, seed=33, mode="indel")],
        "lowcomplex_edge": lambda k: [r[0] for r in make_reads_lowcomplex(st, regions, 20, 100, k, "edge", seed=34)],
        "lowcomplex_inside": lambda k: [r[0] for r in make_reads_lowcomplex(st, regions, 20, 100, k, "inside", seed=35)]}
for name, fn in scen.items():
    for k in (6, 10, 15):
        res = run_case(name, sp, fn(k), k)
        for key, _, desc, _c in ALGOS:
            r = res[key]
            md.append(f"| {name} | {k} | {desc} | {r['hits']:.0f} | {r['windows']:.0f} | {r['bp']:.0f} | {r['ms']:.3f} | {r['sel']:.3f} |")
md.append("\n### Unmodified E. coli 536, uniform edits, m = 100, 100 reads (regression check)\n")
md.append("| k | algorithm | hits | windows | bp verified | query ms | seed selection ms |")
md.append("|---|---|---|---|---|---|---|")
mp = os.path.join(DATA, "ecoli536.txt")
for k in (4, 8, 10, 12, 15):
    reads = [r[0] for r in make_reads(main, 100, 100, k, seed=11)]
    res = run_case("ecoli_uniform", mp, reads, k)
    for key, _, desc, _c in ALGOS:
        r = res[key]
        md.append(f"| {k} | {desc} | {r['hits']:.0f} | {r['windows']:.0f} | {r['bp']:.0f} | {r['ms']:.3f} | {r['sel']:.3f} |")
# before/after summary on the stress case that motivated the change
r = {(x["case"], x["k"], x["algo"]): x for x in rows}
b, a, e = r[("lowcomplex_edge", 10, "filter_fa")], r[("lowcomplex_edge", 10, "filter_fa_w")], r[("lowcomplex_edge", 10, "filter")]
c1k = r[("lowcomplex_edge", 10, "filter_fa_w1k")]
md.append(f"\n**Before/after on the case that verified more bases (lowcomplex_edge, k = 10):** equal split {e['hits']:.0f} hits, {e['windows']:.0f} windows, "
          f"{e['bp']:.0f} bp, {e['ms']:.2f} ms; min-hits DP {b['hits']:.0f} hits, {b['windows']:.0f} windows, {b['bp']:.0f} bp, {b['ms']:.2f} ms; "
          f"hits + lambda * windows DP {a['hits']:.0f} hits, {a['windows']:.0f} windows, {a['bp']:.0f} bp, {a['ms']:.2f} ms (seed selection {a['sel']:.2f} ms); "
          f"same with hit-list cap 1000: {c1k['hits']:.0f} hits, {c1k['windows']:.0f} windows, {c1k['bp']:.0f} bp, {c1k['ms']:.2f} ms (selection {c1k['sel']:.2f} ms).")
for alg in ("filter_fa_w", "filter_fa_w1k"):
    better = sum(1 for (c, k, al), x in r.items() if al == alg and x["ms"] < r[(c, k, "filter_fa")]["ms"] * 0.95)
    worse = sum(1 for (c, k, al), x in r.items() if al == alg and x["ms"] > r[(c, k, "filter_fa")]["ms"] * 1.05)
    tot = sum(1 for (c, k, al) in r if al == alg)
    md.append(f"`{alg}` versus `filter_fa` over all {tot} settings: more than 5% faster in {better}, more than 5% slower in {worse}, within 5% in {tot - better - worse}.")
open(os.path.join(RES, "objective_test.md"), "w").write("\n".join(md) + "\n")
with open(os.path.join(RES, "objective_test.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("\n".join(md[-4:]))
