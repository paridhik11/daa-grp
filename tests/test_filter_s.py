#!/usr/bin/env python3
"""Dedicated stress test for the generalised pigeonhole filters (s = 1, 2, 3; equal-split and frequency-aware seeds).
Reads are planted occurrences with up to k edits (uniform, indel-heavy, clustered), texts include repeats, tandem/low-complexity regions and N,
k from 1 to 25, m from 20 to 200, so true matches exist and the s-hit band logic is exercised (including edits that shift the diagonal by k).
Fallback is disabled so the pure filter path runs.  Also filter_auto with the s-hit options forced on.  Must equal Sellers exactly."""
import os
import random
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from common import run_bin  # noqa: E402
from simulate import inject_low_complexity, make_genome, mutate  # noqa: E402

os.environ["APPROX_NOFALLBACK"] = "1"
rng = random.Random(41)
ALGOS = ["filter_s1", "filter_s2", "filter_s3", "filter_fs2", "filter_fs3"]
bad = cases = planted_found = planted = 0
with tempfile.TemporaryDirectory() as d:
    tp = os.path.join(d, "t.txt")
    for c in range(160):
        n = rng.choice([2000, 8000, 20000])
        g = make_genome(n, seed=rng.randrange(10**6), repeat_frac=rng.choice([0, 0.2]), gc=rng.choice([0.3, 0.5, 0.7]))
        kind = rng.choice(["plain", "lowcomplex", "withN"])
        if kind == "lowcomplex":
            g, _ = inject_low_complexity(g, 0.3, seed=rng.randrange(10**6))
        elif kind == "withN":
            gl = list(g)
            for _ in range(rng.randint(1, n // 200)):
                gl[rng.randrange(n)] = "N"
            g = "".join(gl)
        open(tp, "w").write(g + "\n")
        k = rng.choice([1, 2, 4, 7, 10, 14, 18, 25])
        pats = []
        for _ in range(5):
            m = rng.choice([20, 40, 60, 100, 150, 200])
            m = max(m, k + 3)
            st = rng.randrange(0, n - m + 1)
            mode = rng.choice(["uniform", "indel", "cluster"])
            pats.append(mutate(g[st:st + m].replace("N", "A"), rng.randint(0, k), rng, mode))
        _, ref = run_bin("sellers", tp, pats, k, tag="fs")
        for a in ALGOS:
            _, got = run_bin(a, tp, pats, k, tag="fs")
            cases += 1
            if {i: v for i, v in got.items() if v} != {i: v for i, v in ref.items() if v}:
                bad += 1
                print(f"MISMATCH case={c} algo={a} kind={kind} n={n} k={k} lens={[len(p) for p in pats]}")
        # filter_auto with the s-hit options forced on (trigger 0): whichever procedure it picks must be exact
        for trig in ("0", "30"):
            os.environ["APPROX_AUTOSTRIG"] = trig
            _, got = run_bin("filter_auto", tp, pats, k, tag="fs")
            cases += 1
            if {i: v for i, v in got.items() if v} != {i: v for i, v in ref.items() if v}:
                bad += 1
                print(f"MISMATCH case={c} algo=filter_auto trig={trig} kind={kind} n={n} k={k}")
        os.environ.pop("APPROX_AUTOSTRIG", None)
        planted += sum(1 for i in range(len(pats)) if ref.get(i)); planted_found += sum(1 for i in range(len(pats)) if ref.get(i) and run_bin("filter_s3", tp, pats, k, tag="fs")[1].get(i))
print(f"{cases} comparisons (160 settings, 800 reads, k up to 25, m 20..200, s = 1..3, equal and frequency-aware seeds), {bad} failures; "
      f"reads with a true match: {planted}")
sys.exit(1 if bad else 0)
