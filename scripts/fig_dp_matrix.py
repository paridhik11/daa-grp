#!/usr/bin/env python3
"""Figure: the DP matrix with the Ukkonen band, Landau-Vishkin diagonals, a seed hit and its candidate window.

Small example (m = 24, n = 64, k = 3): P occurs in T with 2 edits.  k + 1 = 4 equal pieces of length 6; one piece occurs
exactly in T (the seed hit).  Drawn: Sellers matrix D (cells with D <= k coloured: the only cells Ukkonen's cutoff must compute),
diagonals d0 - k .. d0 + k around the seed diagonal (the band a window verification looks at), the seed hit on its diagonal,
the candidate window columns [p - a - k, p - a + m + k), and the furthest-reaching points L[d][e] of Landau-Vishkin.
Writes results/dp_matrix.png."""
import os
import random

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
rng = random.Random(4)
m, n, k = 24, 64, 3
T = "".join(rng.choice("ACGT") for _ in range(n))
st = 21
P = list(T[st:st + m])
P[9] = "A" if P[9] != "A" else "C"      # substitution inside piece 1 (touches that piece)
del P[17]                                # deletion inside piece 2
P = "".join(P); m = len(P)
q = m // (k + 1)
pieces = [(i * q, P[i * q:(i + 1) * q]) for i in range(k + 1)]
hits = [(a, [j for j in range(n - q + 1) if T[j:j + q] == pc]) for a, pc in pieces]
D = np.zeros((m + 1, n + 1), dtype=int)
D[:, 0] = np.arange(m + 1)
for j in range(1, n + 1):
    for i in range(1, m + 1):
        D[i, j] = min(D[i - 1, j - 1] + (P[i - 1] != T[j - 1]), D[i - 1, j] + 1, D[i, j - 1] + 1)
a_hit = next((a, p) for a, ps in hits for p in ps if abs(p - a - st) <= k)  # the true-locus seed hit
a0, p0 = a_hit
d0 = p0 - a0
fig, ax = plt.subplots(figsize=(11.5, 7.2))
cmap = plt.get_cmap("Blues_r")
img = np.where(D <= k, D, k + 1).astype(float)
ax.imshow(img, cmap=cmap, vmin=0, vmax=k + 1.6, aspect="equal", origin="upper", extent=(-0.5, n + 0.5, m + 0.5, -0.5))
for i in range(m + 1):
    for j in range(n + 1):
        if D[i, j] <= k:
            ax.text(j, i, str(D[i, j]), ha="center", va="center", fontsize=5.5, color="white" if D[i, j] < 2 else "black")
# candidate window and band
ws, we = max(0, p0 - a0 - k), min(n, p0 - a0 + m + k)
ax.add_patch(Rectangle((ws - 0.5, -0.5), we - ws, m + 1, fill=False, ec="tab:green", lw=2.0, ls="--", label=f"candidate window [{ws}, {we}) = [p - a - k, p - a + m + k)"))
ii = np.arange(0, m + 1)
for dd, ls in ((d0 - k, ":"), (d0 + k, ":"), (d0, "-")):
    ax.plot(ii + dd, ii, color="tab:red", lw=1.2 if dd != d0 else 1.8, ls=ls, label=("band: diagonals d0 - k .. d0 + k" if dd == d0 - k else ("seed diagonal d0 = p - a" if dd == d0 else None)))
ax.plot(np.arange(q) + p0, np.arange(q) + a0, color="gold", lw=5, solid_capstyle="butt", label=f"seed hit: piece at pattern offset a = {a0} found at text position p = {p0}", zorder=5)
# Landau-Vishkin furthest reaching points on the band diagonals (e = k)
NEG = -10 ** 9
L = {}
for e in range(k + 1):
    for d in range(-k, n - m + k + 1):
        if e == 0:
            best = 0 if d >= 0 else NEG
        else:
            best = max(L.get((d, e - 1), NEG) + 1, L.get((d - 1, e - 1), NEG), L.get((d + 1, e - 1), NEG) + 1)
        if best <= NEG // 2 or best + d > n or best + d < 0:
            L[(d, e)] = NEG; continue
        best = min(best, m); i = best
        while i < m and i + d < n and P[i] == T[i + d]:
            i += 1
        L[(d, e)] = i
xs, ys = [], []
for d in range(d0 - k, d0 + k + 1):
    r = L.get((d, k), NEG)
    if r > 0:
        xs.append(r + d); ys.append(r)
ax.scatter(xs, ys, marker="D", s=22, color="k", zorder=6, label="Landau-Vishkin L[d][k]: furthest row reached with k errors")
ax.set_xlabel("text position j (column)"); ax.set_ylabel("pattern position i (row)")
ax.set_title(f"Sellers DP matrix, m = {m}, n = {n}, k = {k}: colored cells have D <= k (the cells Ukkonen's cutoff computes)", fontsize=9)
h_, l_ = ax.get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=2, fontsize=7.5, framealpha=0.95)
ax.set_xlim(-0.5, n + 0.5); ax.set_ylim(m + 0.5, -0.5)
fig.tight_layout(rect=(0, 0.11, 1, 1)); fig.savefig(os.path.join(RES, "dp_matrix.png"), dpi=170)
print("dp_matrix.png written; seed hit a =", a0, "p =", p0, "d0 =", d0, "min D[m][.] =", int(D[m].min()))
