#!/usr/bin/env python3
"""Dedicated stress test for the packed Landau-Vishkin (`lv`) and `lv_direct`: long exact runs (extension across several 32-base blocks),
low-complexity text, tandem repeats, N in text and in the pattern, k from 0 to 40, m up to 300.  Must equal Sellers exactly."""
import os
import random
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from common import run_bin  # noqa: E402
from simulate import inject_low_complexity, make_genome, mutate  # noqa: E402

rng = random.Random(77)
bad = cases = 0
with tempfile.TemporaryDirectory() as d:
    tp = os.path.join(d, "t.txt")
    for c in range(200):
        n = rng.choice([300, 2000, 9000])
        g = make_genome(n, seed=rng.randrange(10**6), repeat_frac=rng.choice([0, 0.3]), gc=rng.choice([0.3, 0.5, 0.7]))
        kind = rng.choice(["plain", "lowcomplex", "tandem", "withN"])
        if kind == "lowcomplex":
            g, _ = inject_low_complexity(g, 0.5, seed=rng.randrange(10**6))
        elif kind == "tandem":
            unit = "".join(rng.choice("ACGT") for _ in range(rng.choice([1, 2, 3, 7])))
            g = (unit * (n // len(unit) + 1))[:n]
            gl = list(g)
            for _ in range(rng.randint(0, 6)):
                gl[rng.randrange(n)] = rng.choice("ACGT")
            g = "".join(gl)
        elif kind == "withN":
            gl = list(g)
            for _ in range(rng.randint(1, max(2, n // 100))):
                gl[rng.randrange(n)] = "N"
            if rng.random() < 0.5:
                gl[0] = "N"; gl[-1] = "N"
            g = "".join(gl)
        open(tp, "w").write(g + "\n")
        k = rng.choice([0, 1, 2, 5, 9, 15, 25, 40])
        pats = []
        for _ in range(4):
            m = rng.choice([1, 5, 31, 32, 33, 64, 100, 150, 300])
            m = min(m, n - 1)
            st = rng.randrange(0, n - m + 1)
            p = mutate(g[st:st + m].replace("N", "A"), rng.randint(0, min(k + 2, m)), rng, rng.choice(["uniform", "indel", "cluster"]))
            if kind == "withN" and rng.random() < 0.5 and len(p) > 2:
                pl = list(p); pl[rng.randrange(len(pl))] = "N"; p = "".join(pl)
            pats.append(p)
        _, ref = run_bin("sellers", tp, pats, k, tag="lvt")
        for a in ("lv", "lv_stream", "lv_direct"):
            _, got = run_bin(a, tp, pats, k, tag="lvt")
            cases += 1
            if {i: s for i, s in got.items() if s} != {i: s for i, s in ref.items() if s}:
                bad += 1
                print(f"MISMATCH case={c} algo={a} kind={kind} n={n} k={k} lens={[len(p) for p in pats]}")
print(f"{cases} comparisons (200 settings, 800 reads, k up to 40, m up to 300, tandem/low-complexity/N texts), {bad} failures")
sys.exit(1 if bad else 0)
