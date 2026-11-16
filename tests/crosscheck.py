#!/usr/bin/env python3
"""Randomised cross-checks (all against Sellers, the DP oracle).

  exact   : banded, myers, lv (Landau-Vishkin), filter, filter_fa must return EXACTLY Sellers' end positions
  heuristic: blast1, blast2 must return a SUBSET
  edlib   : edlib (HW mode, task=locations) reports the ends of the OPTIMAL alignments; it must equal
            our `myers_best` output (ends with minimum distance, if <= k) and that minimum must equal
            Sellers' minimum.  Cases containing N are excluded for edlib (edlib treats N == N).
Scenarios: uniform random, low-complexity text, indel-heavy reads, clustered edits, text/pattern with N.
"""
import os
import random
import subprocess
import sys
import tempfile

import edlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from simulate import make_genome, mutate, inject_low_complexity  # noqa: E402

BIN = os.path.join(ROOT, "approx")


def run(algo, text, pats, k):
    r = subprocess.run([BIN, algo, text, pats, str(k), "--dump"], capture_output=True, text=True, check=True)
    sets = {}
    for line in r.stdout.splitlines():
        if line.startswith("M "):
            _, i, e = line.split()
            sets.setdefault(int(i), set()).add(int(e))
    return sets


def main(cases=300, seed=7):
    rng = random.Random(seed)
    bad = 0
    n_edlib = 0
    stats = {"cases": 0}
    with tempfile.TemporaryDirectory() as d:
        tf, pf = os.path.join(d, "t.txt"), os.path.join(d, "p.txt")
        for c in range(cases):
            n = rng.choice([50, 300, 2000, 8000, 30000])
            g = make_genome(n, seed=rng.randrange(10**6), repeat_frac=rng.choice([0, 0.1]),
                            gc=rng.choice([0.3, 0.5, 0.7]))
            scen = rng.choice(["uniform", "lowcomplex", "indel", "cluster", "withN"])
            if scen == "lowcomplex" and n >= 2000:
                g, _ = inject_low_complexity(g, 0.3, seed=rng.randrange(10**6))
            if scen == "withN":
                gl = list(g)
                for _ in range(max(1, n // 400)):
                    gl[rng.randrange(n)] = "N"
                if rng.random() < 0.5:  # N at the very end (tail handling of the seed index)
                    gl[-1] = "N"
                g = "".join(gl)
            k = rng.choice([0, 1, 2, 3, 5, 8, 13])
            mode = {"indel": "indel", "cluster": "cluster"}.get(scen, "uniform")
            pats = []
            for _ in range(6):
                m = rng.choice([1, 2, 5, 12, 33, 63, 64, 65, 100, 130, 200])
                if rng.random() < 0.7 and m < n:
                    st = rng.randrange(0, n - m + 1)
                    pats.append(mutate(g[st:st + m], rng.randint(0, k + 2), rng, mode))
                else:
                    pats.append("".join(rng.choice("ACGT") for _ in range(m)))
            if scen == "withN":
                pats[0] = pats[0][:len(pats[0]) // 2] + "N" + pats[0][len(pats[0]) // 2 + 1:] if len(pats[0]) > 2 else pats[0]
            open(tf, "w").write(g + "\n")
            open(pf, "w").write("\n".join(pats) + "\n")
            ref = run("sellers", tf, pf, k)
            for algo in ("banded", "myers", "lv", "filter", "filter_fa", "filter_fa_w"):
                got = run(algo, tf, pf, k)
                if got != ref:
                    bad += 1
                    print(f"MISMATCH case={c} scen={scen} algo={algo} n={n} k={k}")
            for algo in ("blast1", "blast2"):
                bl = run(algo, tf, pf, k)
                for i in range(len(pats)):
                    if not bl.get(i, set()) <= ref.get(i, set()):
                        bad += 1
                        print(f"{algo} UNSOUND case={c} scen={scen}")
            if scen != "withN":
                best = run("myers_best", tf, pf, k)
                for i, p in enumerate(pats):
                    r = edlib.align(p, g, mode="HW", task="locations", k=k)
                    ed_ends = {e for _, e in r["locations"]} if r["editDistance"] >= 0 else set()
                    if ed_ends != best.get(i, set()):
                        bad += 1
                        print(f"EDLIB MISMATCH case={c} scen={scen} pat={i} k={k} m={len(p)}")
                    # minimum distance equals Sellers' best: ends in `ref` exist iff edlib found distance <= k
                    if (r["editDistance"] >= 0) != bool(ref.get(i)):
                        bad += 1
                        print(f"EDLIB reachability MISMATCH case={c}")
                    n_edlib += 1
            stats["cases"] += 1
    print(f"{stats['cases']} cases ({n_edlib} edlib pattern checks), {bad} failures")
    return bad


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
