#!/usr/bin/env python3
"""Dedicated stress test for myers_pruned (block pruning): many blocks (m up to 600), k from 0 to well above 64 and above m,
uniform / low-complexity / N-containing text.  Must equal Sellers exactly (also checked against plain Myers)."""
import os
import random
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from common import run_bin  # noqa: E402
from simulate import inject_low_complexity, make_genome, mutate  # noqa: E402

rng = random.Random(33)
bad = cases = 0
with tempfile.TemporaryDirectory() as d:
    tp = os.path.join(d, "t.txt")
    for c in range(160):
        n = rng.choice([1000, 6000, 20000])
        g = make_genome(n, seed=rng.randrange(10**6), repeat_frac=rng.choice([0, 0.1]), gc=rng.choice([0.3, 0.5, 0.7]))
        if rng.random() < 0.3:
            g, _ = inject_low_complexity(g, 0.3, seed=rng.randrange(10**6))
        if rng.random() < 0.25:
            gl = list(g)
            for _ in range(max(1, n // 200)):
                gl[rng.randrange(n)] = "N"
            g = "".join(gl)
        open(tp, "w").write(g + "\n")
        pats = []
        k = rng.choice([0, 1, 5, 30, 63, 64, 65, 100, 129, 150, 200])
        for _ in range(4):
            m = rng.choice([64, 65, 127, 128, 129, 200, 300, 450, 600])
            if m >= n:
                m = n // 2
            if rng.random() < 0.7:
                st = rng.randrange(0, n - m + 1)
                pats.append(mutate(g[st:st + m], rng.randint(0, min(k + 3, m)), rng, rng.choice(["uniform", "indel", "cluster"])))
            else:
                pats.append("".join(rng.choice("ACGT") for _ in range(m)))
        _, ref = run_bin("sellers", tp, pats, k, tag="pr")
        for a in ("myers_pruned", "myers"):
            _, got = run_bin(a, tp, pats, k, tag="pr")
            cases += 1
            if {i: s for i, s in got.items() if s} != {i: s for i, s in ref.items() if s}:
                bad += 1
                print(f"MISMATCH case={c} algo={a} n={n} k={k} lens={[len(p) for p in pats]}")
print(f"{cases} comparisons (160 settings, 640 reads, m up to 600, k up to 200), {bad} failures")
sys.exit(1 if bad else 0)
