#!/usr/bin/env python3
"""Validate alignment output on thousands of random reads.

For every read: (1) mapped iff the Sellers minimum distance is <= k (Sellers at k = d finds a match, at k = d - 1 none);
(2) reported distance == that minimum; (3) CIGAR edit count (X + I + D) == reported distance; (4) the CIGAR replays correctly
against read and text (= means equal and not N, X means different or N, read and text consumption match, within bounds);
(5) end - 1 is one of Sellers' end positions for k.  Writes results/mapping_validation.md and results/mapping_example.txt.
"""
import os
import random
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import common  # noqa: E402
from common import BIN, DATA, RES, run_bin  # noqa: E402
from simulate import make_genome, make_reads, mutate  # noqa: E402


def parse_map(text_path, reads, k):
    pf = os.path.join(DATA, "_map_reads.txt")
    open(pf, "w").write("\n".join(reads) + "\n")
    out = subprocess.run([BIN, "map", text_path, pf, str(k)], capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        if line.startswith("A "):
            t = line.split()
            res[int(t[1])] = t[2:]
    return res


def replay(cigar, read, text, start, end):
    ri, ti, edits = 0, start, 0
    for n, op in re.findall(r"(\d+)([=XID])", cigar):
        n = int(n)
        for _ in range(n):
            if op in "=X":
                if ri >= len(read) or ti >= len(text):
                    return None
                eq = read[ri] == text[ti] and text[ti] != "N"
                if (op == "=") != eq:
                    return None
                if op == "X":
                    edits += 1
                ri += 1; ti += 1
            elif op == "I":
                ri += 1; edits += 1
            else:
                ti += 1; edits += 1
    return (ri == len(read) and ti == end and edits)


def main(runs=320, per=10, seed=21):
    rng = random.Random(seed)
    checked = mapped = unmapped = bad = 0
    with tempfile.TemporaryDirectory() as d:
        tp = os.path.join(d, "t.txt")
        for r in range(runs):
            n = rng.choice([200, 1000, 5000, 20000])
            g = make_genome(n, seed=rng.randrange(10**6), repeat_frac=rng.choice([0, 0.1]), gc=rng.choice([0.3, 0.5, 0.7]))
            if rng.random() < 0.25:
                gl = list(g)
                for _ in range(max(1, n // 300)):
                    gl[rng.randrange(n)] = "N"
                g = "".join(gl)
            open(tp, "w").write(g + "\n")
            k = rng.choice([0, 1, 2, 3, 5, 8, 12])
            mode = rng.choice(["uniform", "indel", "cluster"])
            reads = []
            for _ in range(per):
                m = rng.choice([1, 2, 5, 12, 33, 64, 65, 100, 150])
                if rng.random() < 0.75 and m < n:
                    st = rng.randrange(0, n - m + 1)
                    reads.append(mutate(g[st:st + m], rng.randint(0, k + 2), rng, mode))
                else:
                    reads.append("".join(rng.choice("ACGT") for _ in range(m)))
            reads = [x for x in reads if x]
            res = parse_map(tp, reads, k)
            run_bin("mindist", tp, reads, k, tag="mp")
            dist = dict(common.LAST_DIST)
            _, sel_k = run_bin("sellers", tp, reads, k, tag="mp")
            for i, rd in enumerate(reads):
                checked += 1
                a = res.get(i)
                d_true = dist[i]
                if d_true > k:
                    unmapped += 1
                    if a is None or a[0] != "*":
                        bad += 1; print(f"BAD: expected unmapped run={r} read={i}")
                    continue
                mapped += 1
                if a is None or a[0] == "*":
                    bad += 1; print(f"BAD: expected mapped run={r} read={i}"); continue
                dist_r, start, end, cig = int(a[0]), int(a[1]), int(a[2]), a[3]
                # (1)+(2) distance equals the Sellers minimum: Sellers at k=d finds a match, at k=d-1 none
                rd = reads[i]
                _, s_hi = run_bin("sellers", tp, [rd], dist_r, tag="mp2")
                s_lo = run_bin("sellers", tp, [rd], dist_r - 1, tag="mp2")[1] if dist_r > 0 else {}
                if dist_r != d_true or not s_hi.get(0) or s_lo.get(0):
                    bad += 1; print(f"BAD: distance run={r} read={i} reported={dist_r} min={d_true}")
                ed = replay(cig, rd, g, start, end)
                if ed is None or ed is False or ed != dist_r:
                    # replay returns the edit count (0 is falsy, handle explicitly)
                    if not (ed == 0 and dist_r == 0):
                        bad += 1; print(f"BAD: cigar run={r} read={i} cigar={cig} edits={ed} dist={dist_r}")
                if (end - 1) not in sel_k.get(i, set()):
                    bad += 1; print(f"BAD: end not a Sellers end run={r} read={i}")
    lines = ["## Alignment output validation\n",
             f"{checked} random reads in {runs} runs (texts of 200 to 20,000 bp, some with N; k in 0..12; uniform, indel-heavy and clustered edits; "
             f"{mapped} mapped, {unmapped} unmapped): {bad} failures. Checks per mapped read: distance equals the Sellers minimum (Sellers at k = d matches, at k = d - 1 does not); "
             "CIGAR edit count (X + I + D) equals the reported distance; the CIGAR replays correctly against read and text; the end is a Sellers end position. "
             "Unmapped reads have Sellers minimum distance above k.\n"]
    open(os.path.join(RES, "mapping_validation.md"), "w").write("\n".join(lines) + "\n")
    print(lines[1])
    return bad


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
