#!/usr/bin/env python3
"""Expected-case analysis of the equal-split filter: formula versus measured seed hits (and predicted crossover).

Setting: m = 100, k+1 pieces of length q = floor(m/(k+1)); the seed of each piece is its first s = min(q, 11) characters;
text length n, N = n - s + 1 seed windows, lambda_s = N / 4^s.

(A) Uniform random patterns.  For ANY text, E[hits of one seed] = N / 4^s exactly (every window equals a given uniform
    random s-mer with probability 4^-s), so  E[H] = (k+1) N 4^-s.
(B) Reads sampled from the text with exactly k random edits.  A seed untouched by edits (probability pi ~ (1 - s/m)^k,
    since each edit hits a given s-character seed with probability about s/m) is an exact substring of the text, so its
    expected hits are the size-biased count  C_s = sum_w c(w)^2 / N  (c(w) = occurrences of s-mer w, the seed itself included);
    an edited seed behaves like a random string: N 4^-s.   E[H] ~ (k+1) [ pi C_s + (1 - pi) N 4^-s ].
    On uniform random text C_s = 1 + lambda_s, so this is (k+1)(lambda_s + pi).
    Piece-length detail: edits change the read length to m' = m + (insertions - deletions), so q = floor(m'/(k+1)) and
    s = min(q, 11) differ between reads (a read one base shorter can have a one-shorter seed, which has 4x more hits).
    Each edit is a substitution, insertion or deletion with probability 1/3, so (insertions, deletions) is multinomial and the
    formulas above are averaged over that distribution.
    On a real genome C_s exceeds 1 + lambda_s because of repeats and composition bias: that is the deviation.
Windows and verified bases (hits placed uniformly): lambda_w = H (m + 2k) / n,  W ~ H e^-lambda_w,  V ~ n (1 - e^-lambda_w).
Predicted crossover: smallest k with  c_hit H + c_win W + c_bp B V  >=  a_scan n B  (constants from results/costmodel.json).
Writes results/expected_hits.csv, results/expected_hits.md, results/expected_hits.png.
"""
import csv
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin  # noqa: E402
from simulate import make_genome, make_reads  # noqa: E402

M = 100
KS = [1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 14, 15, 17, 20]
CM = json.load(open(os.path.join(RES, "costmodel.json")))
rnd = os.path.join(DATA, "genome_random.txt")
if not os.path.exists(rnd):
    open(rnd, "w").write(make_genome(4_600_000, seed=7, gc=0.5, repeat_frac=0.0) + "\n")
TEXTS = [("random_text", "Uniform random text (no repeats)", rnd),
         ("synthetic", "Synthetic genome (3% injected repeats)", os.path.join(DATA, "genome_synth.txt")),
         ("ecoli536", "E. coli 536 (real)", os.path.join(DATA, "ecoli536.txt"))]
codes = {"A": 0, "C": 1, "G": 2, "T": 3}


def spectrum_C(arr, s):
    N = len(arr) - s + 1
    code = np.zeros(N, dtype=np.int64)
    for i in range(s):
        code = (code << 2) | arr[i:i + N]
    cnt = np.bincount(code, minlength=4 ** s).astype(np.float64)
    return float((cnt ** 2).sum() / N)


def exp_hits(k, n, Cs, spectrum):
    """E[hits] for reads with exactly k edits (each sub/ins/del with probability 1/3), averaged over read length m' = m + ins - del."""
    tot = 0.0
    f = math.factorial
    for i in range(k + 1):
        for d in range(k + 1 - i):
            pr = f(k) / (f(i) * f(d) * f(k - i - d)) / 3 ** k
            mlen = M + i - d
            q = mlen // (k + 1)
            if q < 1:
                continue
            s = min(q, 11)
            lam = (n - s + 1) / 4 ** s
            pi = (1 - s / mlen) ** k
            tot += pr * (k + 1) * ((pi * Cs[s] + (1 - pi) * lam) if spectrum else (lam + pi))
    return tot


def crossover(ks, pred_f, pred_s):
    for (k0, a0, b0), (k1, a1, b1) in zip(zip(ks, pred_f, pred_s), list(zip(ks, pred_f, pred_s))[1:]):
        if a0 < b0 <= a1 or (a0 < b0 and a1 >= b1):
            r0, r1 = math.log(a0 / b0), math.log(a1 / b1)
            return k0 + (k1 - k0) * (0 - r0) / (r1 - r0)
    return None


rows, md = [], ["## Expected-case analysis: formula versus measured seed hits (equal split, m = 100)\n"]
md.append("`hits` = seed hits per read, mean over reads. (A) random uniform patterns: formula (k+1) N 4^-s. (B) reads sampled from the text with exactly k edits: "
          "uniform-text formula (k+1)(lambda_s + pi) and spectrum-corrected formula (k+1)[pi C_s + (1-pi) N 4^-s], pi = (1 - s/m')^k, both averaged over the read-length distribution m' = m + insertions - deletions (piece-length detail).\n")
fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
summary = []
for ti, (key, label, path) in enumerate(TEXTS):
    if not os.path.exists(path):
        continue
    txt = open(path).read().strip()
    n = len(txt)
    arr = np.frombuffer(txt.encode(), dtype=np.uint8)
    lut = np.zeros(256, dtype=np.int64)
    for c, v in codes.items():
        lut[ord(c)] = v
    arr = lut[arr]
    rng = np.random.default_rng(11)
    Cs = {}
    md.append(f"\n### {label} (n = {n:,})\n")
    md.append("| k | q | s | (A) formula | (A) measured | (A) ratio | (B) uniform formula | (B) spectrum-corrected | (B) measured | measured / uniform | measured / corrected | C_s | 1 + lambda_s |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    meas_A, form_A, meas_B, form_B, form_Bc = [], [], [], [], []
    for k in KS:
        q = M // (k + 1)
        s = min(q, 11)
        N = n - s + 1
        lam = N / 4 ** s
        for s_ in range(1, 12):
            if s_ not in Cs:
                Cs[s_] = spectrum_C(arr, s_)
        fA = (k + 1) * lam
        fB = exp_hits(k, n, Cs, spectrum=False)
        fBc = exp_hits(k, n, Cs, spectrum=True)
        nreads = 300 if k <= 12 else 60
        rpat = ["".join("ACGT"[v] for v in rng.integers(0, 4, M)) for _ in range(nreads)]
        info, _ = run_bin("filter", path, rpat, k, nofallback=True, dump=False, tag="eh")
        mA = info["hits"] / nreads
        reads = [r[0] for r in make_reads(txt, nreads, M, k, seed=500 + k)]
        info, _ = run_bin("filter", path, reads, k, nofallback=True, dump=False, tag="eh")
        mB = info["hits"] / nreads
        rows.append(dict(text=key, k=k, q=q, s=s, formula_A=fA, measured_A=mA, formula_B=fB, formula_B_corrected=fBc, measured_B=mB, C_s=Cs[s], lam=lam))
        md.append(f"| {k} | {q} | {s} | {fA:.1f} | {mA:.1f} | {mA / fA:.2f} | {fB:.1f} | {fBc:.1f} | {mB:.1f} | {mB / fB:.2f} | {mB / fBc:.2f} | {Cs[s]:.2f} | {1 + lam:.2f} |")
        meas_A.append(mA); form_A.append(fA); meas_B.append(mB); form_B.append(fB); form_Bc.append(fBc)
    ratio = lambda a, b: sorted(abs(math.log(x / y)) for x, y in zip(a, b))
    rA, rB, rBc = ratio(meas_A, form_A), ratio(meas_B, form_B), ratio(meas_B, form_Bc)
    med = lambda r: math.exp(r[len(r) // 2])
    summary.append((label, med(rA), math.exp(rA[-1]), med(rB), math.exp(rB[-1]), med(rBc), math.exp(rBc[-1])))
    # predicted crossover from the formula (uniform-text H) and from the spectrum-corrected H
    B = (M + 63) // 64
    kk = list(range(1, 41))
    def pred_time(Hf, k):
        lw = Hf * (M + 2 * k) / n
        W = Hf * math.exp(-lw); V = n * (1 - math.exp(-lw))
        return (CM["c_hit"] * Hf + CM["c_win"] * W + CM["c_bp"] * B * V) / 1e6
    pf_u, pf_c = [], []
    for k in kk:
        for s_ in range(1, 12):
            if s_ not in Cs:
                Cs[s_] = spectrum_C(arr, s_)
        pf_u.append(pred_time(exp_hits(k, n, Cs, False), k)); pf_c.append(pred_time(exp_hits(k, n, Cs, True), k))
    pscan = [CM["a_scan"] * n * B / 1e6] * len(kk)
    cu, cc = crossover(kk, pf_u, pscan), crossover(kk, pf_c, pscan)
    md.append(f"\nPredicted crossover with a Myers scan (cost-model constants, formula hits): uniform-text formula k = {cu:.1f} (k/m = {cu / 100:.3f})" if cu else "\nPredicted crossover (uniform-text formula): none up to k = 40")
    md[-1] += f"; spectrum-corrected formula k = {cc:.1f} (k/m = {cc / 100:.3f})." if cc else "; spectrum-corrected: none up to k = 40."
    ax = axs[ti]
    ax.plot(KS, meas_A, "o", color="tab:blue", label="measured, random patterns")
    ax.plot(KS, form_A, "-", color="tab:blue", label="formula (A)")
    ax.plot(KS, meas_B, "s", color="tab:red", label="measured, reads from text")
    ax.plot(KS, form_B, "--", color="tab:red", label="uniform-text formula (B)")
    ax.plot(KS, form_Bc, ":", color="tab:red", lw=2, label="spectrum-corrected (B)")
    ax.set_yscale("log"); ax.set_xlabel("k (m = 100)"); ax.set_ylabel("seed hits per read"); ax.set_title(label, fontsize=9)
    ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7)
fig.suptitle("Equal-split filter: expected versus measured seed hits")
fig.tight_layout(); fig.savefig(os.path.join(RES, "expected_hits.png"), dpi=160); plt.close(fig)
md.append("\n### Agreement (factor error over k = " + ", ".join(map(str, KS)) + ": median / worst)\n")
md.append("| text | (A) random patterns vs formula | (B) reads vs uniform formula | (B) reads vs spectrum-corrected formula |")
md.append("|---|---|---|---|")
for lab, a, aw, b, bw, c, cw in summary:
    md.append(f"| {lab} | {a:.2f}x / {aw:.2f}x | {b:.2f}x / {bw:.2f}x | {c:.2f}x / {cw:.2f}x |")
with open(os.path.join(RES, "expected_hits.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
open(os.path.join(RES, "expected_hits.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md[-8:]))
