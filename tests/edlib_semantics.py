#!/usr/bin/env python3
"""Shows, on a tiny case, what Edlib HW mode returns versus the all-ends-within-k definition used by Sellers.

Edlib returns the MINIMUM distance and only the end positions achieving it. Our k-differences definition
reports every end with distance <= k.  So Edlib is compared with `myers_best` (min-distance ends), not
with the full Sellers set.
"""
import os
import subprocess
import sys
import tempfile

import edlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "approx")
T, P, K = "TTACGTTTTACGATT", "ACGT", 1


def ends(algo):
    with tempfile.TemporaryDirectory() as d:
        open(d + "/t", "w").write(T + "\n"); open(d + "/p", "w").write(P + "\n")
        out = subprocess.run([BIN, algo, d + "/t", d + "/p", str(K), "--dump"], capture_output=True, text=True).stdout
    return sorted(int(l.split()[2]) for l in out.splitlines() if l.startswith("M "))


r = edlib.align(P, T, mode="HW", task="locations", k=K)
print("text   :", T)
print("pattern:", P, " k =", K)
print("Sellers  (all ends with distance <= k):", ends("sellers"))
print("myers_best (min-distance ends only)   :", ends("myers_best"))
print("Edlib HW  distance:", r["editDistance"], " ends:", sorted(e for _, e in r["locations"]))
assert sorted(e for _, e in r["locations"]) == ends("myers_best")
assert len(ends("sellers")) > len(ends("myers_best")), "demo case must show a difference"
print("OK: Edlib == min-distance ends, which is a strict subset of the all-within-k set here.")
