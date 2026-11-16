#!/usr/bin/env python3
"""All benchmarks -> results/bench.csv.   Usage: python3 scripts/run_bench.py [--quick]

Experiments (main text = data/genome.txt, the largest real bacterial genome found):
  error    : m=100, k sweep                      (all algorithms)
  length   : m sweep at 5% error
  textsize : prefixes of the main text, m=100, k=5
  datasets : every real genome, m=100, k=5
  stress   : uniform / clustered / indel-heavy / low-complexity-inside / low-complexity-edge reads
Ground truth = Myers (identical to Sellers: asserted on the patterns Sellers runs, checked in tests/crosscheck.py).
Every exact algorithm's match set is asserted equal to the truth; heuristics are asserted to be subsets.
Edlib and WFA2 return the BEST match per read (minimum distance, and only the ends at that minimum), not all
ends within k.  Rows carry semantics="best-per-read" and are NOT comparable to all-matches-within-k rows on
output, only on time.  Edlib is checked for equality with our `myers_best` (min-distance ends); WFA2 for
validity (its end is a minimum-distance end).  Timings include Python-binding call overhead (~0.1 ms).
"""
import csv
import json
import os
import sys
import time

import numpy as np

import edlib
import pywfa

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin  # noqa: E402
from simulate import (inject_low_complexity, make_reads, make_reads_lowcomplex)  # noqa: E402

QUICK = "--quick" in sys.argv
NP = {"sellers": 1, "banded": 20, "myers": 20, "edlib": 20, "wfa": 5, "filter": 20, "filter_nofb": 20,
      "filter_fa": 20, "filter_fa_nofb": 20, "blast1": 20, "blast2": 20}
if QUICK:
    NP = {a: min(v, 3) for a, v in NP.items()}
BIN_ALGOS = {"sellers": ("sellers", False), "banded": ("banded", False), "myers": ("myers", False),
             "filter": ("filter", False), "filter_nofb": ("filter", True),
             "filter_fa": ("filter_fa", False), "filter_fa_nofb": ("filter_fa", True),
             "blast1": ("blast1", False), "blast2": ("blast2", False)}
EXACT = {"sellers", "banded", "myers", "filter", "filter_nofb", "filter_fa", "filter_fa_nofb"}
HEUR = {"blast1", "blast2"}
ROWS = []
TEXT_CACHE = {}


def text_of(name):
    if name not in TEXT_CACHE:
        TEXT_CACHE[name] = open(os.path.join(DATA, name + ".txt")).read().strip()
    return TEXT_CACHE[name]


def prefix_file(name, n):
    full = text_of(name)
    if n >= len(full):
        return os.path.join(DATA, name + ".txt"), full
    p = os.path.join(DATA, f"_prefix_{name}_{n}.txt")
    open(p, "w").write(full[:n] + "\n")
    return p, full[:n]


def measure(exp, x, dataset, text_path, text, reads, k, algos, mode="uniform", extra=None):
    n, m = len(text), len(reads[0])
    tb = text.encode()
    truth = None
    results = {}
    # truth first
    info_t, truth = run_bin("myers", text_path, reads[:NP["myers"]], k, tag="t")
    ntruth = {i: truth.get(i, set()) for i in range(min(len(reads), NP["myers"]))}
    for algo in algos:
        npat = min(NP[algo], len(reads))
        pats = reads[:npat]
        row = dict(exp=exp, x=x, dataset=dataset, n=n, m=m, k=k, mode=mode, algo=algo, npat=npat,
                   index_ms=0.0, hits_pp=0, windows_pp=0, verified_bp_pp=0, fallbacks=0, filtered=0,
                   pred_filter_ms_pp=0, pred_scan_ms_pp=0, gen_ms_pp=0, ver_ms_pp=0, sel_ms_pp=0,
                   matches=0, recall_end="", recall_read="", exact_ok="", cells=0, cells_per_s=0,
                   semantics="best-per-read" if algo in ("edlib", "wfa") else "all-within-k")
        if algo in BIN_ALGOS:
            ba, nofb = BIN_ALGOS[algo]
            info, sets = run_bin(ba, text_path, pats, k, nofallback=nofb, tag="b")
            t_ms = info["time_ms"]
            row.update(index_ms=info["index_ms"], hits_pp=info["hits"] / npat, windows_pp=info["windows"] / npat,
                       verified_bp_pp=info["verified_bp"] / npat, fallbacks=int(info["fallbacks"]),
                       filtered=int(info["filtered"]), pred_filter_ms_pp=info["pred_filter_ms"] / npat,
                       pred_scan_ms_pp=info["pred_scan_ms"] / npat, gen_ms_pp=info["gen_ms"] / npat,
                       ver_ms_pp=info["ver_ms"] / npat, sel_ms_pp=info["sel_ms"] / npat,
                       matches=int(info["matches"]))
            if algo == "banded":
                row["cells"] = info["cells_banded"]
            elif algo in ("sellers", "myers"):
                row["cells"] = m * n * npat  # DP cells represented (Myers: 64 per word operation)
            if row["cells"]:
                row["cells_per_s"] = row["cells"] / (t_ms / 1e3)
            if algo in EXACT:
                ref = {i: ntruth.get(i, set()) for i in range(npat)}
                got = {i: sets.get(i, set()) for i in range(npat)}
                row["exact_ok"] = int(got == ref)
                if not row["exact_ok"]:
                    raise SystemExit(f"EXACT MISMATCH {exp} x={x} k={k} {algo}")
            else:
                tot = sum(len(ntruth.get(i, ())) for i in range(npat))
                hit = sum(len(sets.get(i, set()) & ntruth.get(i, set())) for i in range(npat))
                assert all(sets.get(i, set()) <= ntruth.get(i, set()) for i in range(npat)), f"{algo} unsound"
                reads_t = [i for i in range(npat) if ntruth.get(i)]
                row["recall_end"] = hit / tot if tot else 1.0
                row["recall_read"] = (sum(1 for i in reads_t if sets.get(i)) / len(reads_t)) if reads_t else 1.0
        elif algo == "edlib":
            t0 = time.perf_counter(); out = []
            for p in pats:
                out.append(edlib.align(p, tb, mode="HW", task="locations", k=k))
            t_ms = (time.perf_counter() - t0) * 1e3
            _, best = run_bin("myers_best", text_path, pats, k, tag="e")
            for i, r in enumerate(out):
                ends = {e for _, e in r["locations"]} if r["editDistance"] >= 0 else set()
                tr = ntruth.get(i, set())
                assert (r["editDistance"] >= 0) == bool(tr) and ends == best.get(i, set()) and ends <= tr, \
                    f"edlib != min-distance ends {exp} {x}"
            row["matches"] = sum(len(r["locations"]) for r in out if r["editDistance"] >= 0)
            row["exact_ok"] = 1
        elif algo == "wfa":
            t0 = time.perf_counter(); out = []
            for p in pats:
                a = pywfa.WavefrontAligner(p, span="ends-free", pattern_begin_free=0, pattern_end_free=0,
                                           text_begin_free=n, text_end_free=n, distance="levenshtein")
                out.append(a(text))
            t_ms = (time.perf_counter() - t0) * 1e3
            for i, r in enumerate(out):
                tr = ntruth.get(i, set())
                if r.score <= k:
                    assert (r.text_end - 1) in tr, f"wfa end not a valid match {exp} {x}"
                else:
                    assert not tr, f"wfa missed within-k match {exp} {x}"
            row["matches"] = len(out)
            row["exact_ok"] = 1
        row["ms_per_pattern"] = t_ms / npat
        if extra:
            row.update(extra)
        ROWS.append(row)
        print(f"[{exp}] x={x} {dataset} k={k} {algo:15s} {row['ms_per_pattern']:10.3f} ms/pat "
              f"hits/pat={row['hits_pp']:.0f} fb={row['fallbacks']} rec_read={row['recall_read']} ", flush=True)


ALL = ["sellers", "banded", "myers", "edlib", "wfa", "filter", "filter_nofb", "filter_fa", "filter_fa_nofb", "blast1", "blast2"]


HEADLINE_ALGOS = ["myers", "filter", "filter_nofb", "filter_fa", "filter_fa_nofb"]
HEADLINE_FIELDS = ["k", "algo", "npat", "index_ms", "t_med", "t_p10", "t_p90", "t_mean", "hits_med", "hits_p10", "hits_p90",
                   "bp_med", "bp_p10", "bp_p90", "predf_med", "preds_med", "fallbacks", "exact_ok"]


def headline(text_path, text, ks, nreads):
    """E. coli 536, m = 100: filter vs filter_fa vs Myers with `nreads` reads per setting.
    Per-read time / hits / verified bp -> median and 10th/90th percentiles (results/headline.csv, ms)."""
    out = []
    for k in ks:
        reads = [r[0] for r in make_reads(text, nreads, 100, k, seed=11)]
        info_t, truth, tq_t = run_bin("myers", text_path, reads, k, tag="h", times=True)
        for algo in HEADLINE_ALGOS:
            ba, nofb = {"myers": ("myers", False), "filter": ("filter", False), "filter_nofb": ("filter", True),
                        "filter_fa": ("filter_fa", False), "filter_fa_nofb": ("filter_fa", True)}[algo]
            info, sets, tq = run_bin(ba, text_path, reads, k, nofallback=nofb, tag="h", times=True)
            ok = all(sets.get(i, set()) == truth.get(i, set()) for i in range(len(reads)))
            if not ok:
                raise SystemExit(f"HEADLINE EXACT MISMATCH k={k} {algo}")
            ms = np.array([q["ns"] for q in tq]) / 1e6
            h = np.array([q["hits"] for q in tq]); bp = np.array([q["bp"] for q in tq])
            pf = np.array([q["pf"] for q in tq]) / 1e6; ps = np.array([q["ps"] for q in tq]) / 1e6
            P = lambda a, q: float(np.percentile(a, q))
            out.append(dict(k=k, algo=algo, npat=len(reads), index_ms=info["index_ms"], t_med=P(ms, 50), t_p10=P(ms, 10),
                            t_p90=P(ms, 90), t_mean=float(ms.mean()), hits_med=P(h, 50), hits_p10=P(h, 10), hits_p90=P(h, 90),
                            bp_med=P(bp, 50), bp_p10=P(bp, 10), bp_p90=P(bp, 90), predf_med=P(pf, 50), preds_med=P(ps, 50),
                            fallbacks=int(info["fallbacks"]), exact_ok=1))
            print(f"[headline] k={k} {algo:15s} median {P(ms, 50):9.3f} ms  p10-p90 [{P(ms, 10):.3f}, {P(ms, 90):.3f}]  "
                  f"hits(med)={P(h, 50):.0f} fb={int(info['fallbacks'])}", flush=True)
    with open(os.path.join(RES, "headline.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADLINE_FIELDS)
        w.writeheader(); w.writerows(out)


def main():
    main_name = "ecoli536" if os.path.exists(os.path.join(DATA, "ecoli536.txt")) else "genome_synth"
    src = open(os.path.join(DATA, "SOURCES.txt")).read().split("\n")[0]
    if "SYNTHETIC" in src:
        main_name = "genome_synth"
    print("main text:", main_name, flush=True)
    N = len(text_of(main_name))
    tp, tx = prefix_file(main_name, N)

    ks = [1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 14, 15, 17, 20, 30] if not QUICK else [2, 8, 15]
    for k in ks:
        reads = [r[0] for r in make_reads(tx, 20, 100, k, seed=11)]
        measure("error", k / 100, main_name, tp, tx, reads, k, ALL)
    headline(tp, tx, ks, 100 if not QUICK else 10)
    for m in ((50, 100, 200, 400, 800) if not QUICK else (50, 200)):
        k = max(1, m // 20)
        reads = [r[0] for r in make_reads(tx, 20, m, k, seed=12)]
        measure("length", m, main_name, tp, tx, reads, k,
                ["sellers", "banded", "myers", "edlib", "wfa", "filter", "filter_fa", "blast1", "blast2"])
    for n in ([250_000, 500_000, 1_000_000, 2_000_000, N] if not QUICK else [250_000, 1_000_000]):
        p2, t2 = prefix_file(main_name, n)
        reads = [r[0] for r in make_reads(t2, 20, 100, 5, seed=13)]
        measure("textsize", n, main_name, p2, t2, reads, 5,
                ["sellers", "banded", "myers", "edlib", "wfa", "filter", "filter_fa", "blast1", "blast2"])
    for ds in ("ecoli536", "sarscov2", "lambda", "mt_human") if main_name != "genome_synth" else ("genome_synth",):
        if not os.path.exists(os.path.join(DATA, ds + ".txt")):
            continue
        p2, t2 = prefix_file(ds, len(text_of(ds)))
        for k in (5, 10):
            reads = [r[0] for r in make_reads(t2, 20, 100, k, seed=14)]
            measure("datasets", ds, ds, p2, t2, reads, k,
                    ["sellers", "myers", "edlib", "wfa", "filter", "filter_fa", "blast1", "blast2"])
    # ---- stress: main text with ~1% microsatellite / homopolymer arrays injected
    st_text, regions = inject_low_complexity(tx, 0.01, seed=5)
    sp = os.path.join(DATA, "stress_text.txt")
    open(sp, "w").write(st_text + "\n")
    json.dump(regions, open(os.path.join(DATA, "stress_regions.json"), "w"))
    scen = {
        "uniform": lambda k: [r[0] for r in make_reads(st_text, 20, 100, k, seed=31)],
        "clustered": lambda k: [r[0] for r in make_reads(st_text, 20, 100, k, seed=32, mode="cluster")],
        "indel_heavy": lambda k: [r[0] for r in make_reads(st_text, 20, 100, k, seed=33, mode="indel")],
        "lowcomplex_edge": lambda k: [r[0] for r in make_reads_lowcomplex(st_text, regions, 20, 100, k, "edge", seed=34)],
        "lowcomplex_inside": lambda k: [r[0] for r in make_reads_lowcomplex(st_text, regions, 20, 100, k, "inside", seed=35)],
    }
    for name, fn in scen.items():
        for k in ((6, 10, 15) if not QUICK else (6,)):
            measure("stress", name, main_name + "+lowcomplex", sp, st_text, fn(k), k,
                    ["myers", "edlib", "filter", "filter_nofb", "filter_fa", "filter_fa_nofb", "blast1", "blast2"],
                    mode=name)
    fields = list(ROWS[0].keys())
    for r in ROWS:
        for f in r:
            if f not in fields:
                fields.append(f)
    with open(os.path.join(RES, "bench.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(ROWS)
    print("wrote results/bench.csv with", len(ROWS), "rows")


if __name__ == "__main__":
    main()
