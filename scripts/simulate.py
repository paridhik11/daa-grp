#!/usr/bin/env python3
"""Synthetic genome + read simulator (real NCBI genomes are unreachable from the sandbox).

genome: random ACGT (GC-biased) with injected diverged repeats, so filters face realistic
        repeat content.  reads: random substrings with EXACTLY `edits` edits (sub/ins/del),
        so the true edit distance is <= edits and the true location is known.
"""
import argparse
import random

import numpy as np

ALPHA = "ACGT"


def make_genome(n, seed=1, gc=0.5, repeat_frac=0.03, div=0.03):
    rng = np.random.default_rng(seed)
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    g = rng.choice(4, size=n, p=p).astype(np.uint8)
    r = random.Random(seed)
    injected = 0
    while injected < repeat_frac * n:
        L = r.randint(300, 5000)
        if L >= n // 4:
            break
        src = r.randrange(0, n - L)
        dst = r.randrange(0, n - L)
        seg = g[src:src + L].copy()
        for i in range(L):  # diverge by substitutions only (keeps length)
            if r.random() < div:
                seg[i] = (seg[i] + r.randint(1, 3)) % 4
        g[dst:dst + L] = seg
        injected += L
    return "".join(ALPHA[x] for x in g)


def _edit_kind(mode, rng):
    if mode == "indel":  # indel-heavy: 10% substitutions, 90% indels
        return rng.choices("sid", weights=[10, 45, 45])[0]
    return rng.choice("sid")


def mutate(s, edits, rng, mode="uniform"):
    """Apply exactly `edits` random edits (sub / ins / del).
    mode: uniform (equal kinds, random places) | indel (90% indels) | cluster (all edits inside a short window)."""
    s = list(s)
    lo, hi = 0, len(s)
    if mode == "cluster":
        width = min(len(s), max(2 * edits, 12))
        lo = rng.randrange(0, len(s) - width + 1)
        hi = lo + width
    for _ in range(edits):
        kind = _edit_kind(mode, rng)
        hi = min(hi, len(s))
        i = rng.randrange(lo, max(lo + 1, hi))
        if kind == "s":
            s[i] = rng.choice([c for c in ALPHA if c != s[i]])
        elif kind == "i":
            s.insert(i, rng.choice(ALPHA))
        elif len(s) > 1:
            del s[i]
    return "".join(s)


UNITS = ["CA", "CAG", "AT", "A", "GAA", "CCG", "TTAGGG", "GT", "C"]


def inject_low_complexity(genome, frac=0.01, seed=3):
    """Overwrite ~frac of the genome with microsatellite / homopolymer arrays (length preserved).
    Returns (new_genome, [(start, end, unit)])."""
    r = random.Random(seed)
    g = list(genome)
    regions, used, tries = [], 0, 0
    while used < frac * len(g) and tries < 10000:
        tries += 1
        L = r.randint(min(300, len(g) // 8), min(3000, max(1, len(g) // 4)))
        st = r.randrange(0, len(g) - L)
        if any(not (st + L < a or st > b) for a, b, _ in regions):
            continue
        u = r.choice(UNITS)
        g[st:st + L] = list((u * (L // len(u) + 1))[:L])
        regions.append((st, st + L, u))
        used += L
    return "".join(g), sorted(regions)


def make_reads_lowcomplex(genome, regions, count, m, edits, where="inside", seed=2):
    """Reads that lie inside a low-complexity array ('inside') or straddle its boundary ('edge')."""
    rng = random.Random(seed)
    out = []
    for _ in range(count):
        a, b, _u = rng.choice([x for x in regions if x[1] - x[0] >= 2 * m] or regions)
        if where == "inside":
            st = rng.randrange(a, max(a + 1, b - m))
        else:  # half the read in the array, half in unique sequence
            st = (a - m // 2) if rng.random() < 0.5 else (b - m // 2)
        st = max(0, min(len(genome) - m, st))
        out.append((mutate(genome[st:st + m], edits, rng), st, st + m))
    return out


def make_reads(genome, count, m, edits, seed=2, mode="uniform"):
    """Returns list of (read, true_start, true_end_exclusive)."""
    rng = random.Random(seed)
    out = []
    for _ in range(count):
        st = rng.randrange(0, len(genome) - m)
        out.append((mutate(genome[st:st + m], edits, rng, mode), st, st + m))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=4_600_000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="data/genome_synth.txt")
    a = ap.parse_args()
    g = make_genome(a.n, a.seed)
    open(a.out, "w").write(g + "\n")
    print(f"wrote {a.out}: {len(g)} bp")
