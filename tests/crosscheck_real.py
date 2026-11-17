#!/usr/bin/env python3
"""Cross-check every exact algorithm against Sellers on REAL genomes (prefixes of E. coli 536 and B. fragilis, whole
SARS-CoV-2, lambda, human mtDNA) with simulated reads of several kinds, and also mt_orang segments against mt_human
(real homologous sequences).  Exact methods must equal Sellers exactly; heuristics must be a subset."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from common import DATA, run_bin  # noqa: E402
from simulate import make_reads  # noqa: E402

EXACT = ["banded", "myers", "myers_pruned", "lv", "filter", "filter_fa", "filter_fa_w", "filter_auto", "filter_fa_full"]
HEUR = ["blast1", "blast2"]
bad = 0
checks = 0


def text(name, n=None):
    t = open(os.path.join(DATA, name + ".txt")).read().strip()
    return t[:n] if n else t


def check(label, tpath, reads, k):
    global bad, checks
    _, ref = run_bin("sellers", tpath, reads, k, tag="cr")
    for a in EXACT:
        for nofb in ((False, True) if a.startswith("filter") else (False,)):
            _, got = run_bin(a, tpath, reads, k, nofallback=nofb, tag="cr")
            checks += 1
            if got != {i: ref.get(i, set()) for i in range(len(reads)) if ref.get(i)} and got != ref:
                if {i: s for i, s in got.items() if s} != {i: s for i, s in ref.items() if s}:
                    bad += 1
                    print(f"MISMATCH {label} k={k} {a} nofb={nofb}")
    for a in HEUR:
        _, got = run_bin(a, tpath, reads, k, tag="cr")
        checks += 1
        if not all(got.get(i, set()) <= ref.get(i, set()) for i in range(len(reads))):
            bad += 1
            print(f"UNSOUND {label} k={k} {a}")


for name, n in (("ecoli536", 300_000), ("bfragilis", 300_000), ("sarscov2", None), ("lambda", None), ("mt_human", None)):
    p = os.path.join(DATA, name + ".txt")
    if not os.path.exists(p):
        continue
    t = text(name, n)
    tp = os.path.join(DATA, f"_cr_{name}.txt")
    open(tp, "w").write(t + "\n")
    for mode in ("uniform", "cluster", "indel"):
        for k in (3, 8, 15):
            reads = [r[0] for r in make_reads(t, 6, 100, k, seed=hash((name, mode, k)) % 10_000, mode=mode)]
            check(f"{name}/{mode}", tp, reads, k)
# real homologous sequences: random segments of mt_orang against mt_human
if os.path.exists(os.path.join(DATA, "mt_orang.txt")):
    import random
    h, o = text("mt_human"), text("mt_orang")
    tp = os.path.join(DATA, "_cr_mt_human.txt")
    open(tp, "w").write(h + "\n")
    rng = random.Random(5)
    segs = [o[s:s + 100] for s in (rng.randrange(0, len(o) - 100) for _ in range(8))]
    for k in (10, 20, 30):
        check("mt_orang->mt_human", tp, segs, k)
print(f"{checks} comparisons on real data, {bad} failures")
sys.exit(1 if bad else 0)
