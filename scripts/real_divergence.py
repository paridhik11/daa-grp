#!/usr/bin/env python3
"""Real homologous sequences: random segments of mt_orang matched against mt_human, and the reverse.

These are REAL homologs (human versus orangutan mitochondrial DNA), not simulated edits: the true edit distance of a
segment to its best match is whatever evolution produced.  For each segment we compute the minimum edit distance (Myers,
verified against Sellers: Sellers at k = d finds a match and at k = d - 1 finds none), then for several k we measure
recall of the heuristics, time of every method and the filter's candidate counts.
Writes results/divergence.md, results/divergence.png and results/divergence.csv.
"""
import csv
import os
import random
import sys
import time

import edlib
import numpy as np
import pywfa
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
from common import DATA, RES, run_bin  # noqa: E402

NSEG = 200
KS = [10, 15, 20, 25, 30, 40]
ALGOS = ["sellers", "lv", "banded", "myers", "edlib", "wfa", "filter", "filter_nofb", "filter_fa", "filter_fa_nofb", "blast1", "blast2"]
EXACT = {"sellers", "lv", "banded", "myers", "filter", "filter_nofb", "filter_fa", "filter_fa_nofb"}


def seq(name):
    return open(os.path.join(DATA, name + ".txt")).read().strip()


human, orang = seq("mt_human"), seq("mt_orang")
md = ["## Real divergence: mt_orang segments against mt_human, and the reverse\n",
      "Both are real homologous mitochondrial genomes (human, 16,569 bp; orangutan, 16,499 bp), so the edits are real evolutionary differences, "
      "not simulated ones. Segments are random substrings of the query genome; the text is the other genome (forward strand). "
      "The true edit distance is the minimum over all text positions, computed with Myers and verified with Sellers.\n"]
rows = []
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
for di, (qname, tname, q, t) in enumerate((("mt_orang", "mt_human", orang, human), ("mt_human", "mt_orang", human, orang))):
    tpath = os.path.join(DATA, tname + ".txt")
    for m in (100, 200):
        rng = random.Random(100 + m + di)
        segs = [q[s:s + m] for s in (rng.randrange(0, len(q) - m) for _ in range(NSEG))]
        run_bin("mindist", tpath, segs, m, tag="dv")
        dist = dict(common.LAST_DIST)
        d = np.array([dist[i] for i in range(NSEG)])
        # verify with Sellers: match at k = d, none at k = d - 1
        bad = 0
        for dv in sorted(set(d.tolist())):
            idx = [i for i in range(NSEG) if d[i] == dv]
            _, s_hi = run_bin("sellers", tpath, [segs[i] for i in idx], dv, tag="dv")
            _, s_lo = run_bin("sellers", tpath, [segs[i] for i in idx], dv - 1, tag="dv") if dv > 0 else (None, {})
            bad += sum(1 for j in range(len(idx)) if not s_hi.get(j)) + sum(1 for j in range(len(idx)) if s_lo.get(j))
        pct = np.percentile(d / m, [10, 50, 90])
        md.append(f"\n### {qname} segments -> {tname} text, m = {m} ({NSEG} segments)\n")
        md.append(f"True best edit distance d* (Sellers-verified, {bad} verification failures): median {np.median(d):.0f} "
                  f"({100 * pct[1]:.1f}% of m), 10th to 90th percentile {np.percentile(d, 10):.0f} to {np.percentile(d, 90):.0f} "
                  f"({100 * pct[0]:.1f}% to {100 * pct[2]:.1f}%), min {d.min()}, max {d.max()}; "
                  f"{int((d == 0).sum())} segments match exactly.\n")
        if m == 100:
            axs[di][0].hist(d, bins=range(0, int(d.max()) + 2), color="tab:blue" if di == 0 else "tab:orange")
            axs[di][0].set_title(f"{qname} (m=100) vs {tname}: true edit distance (Sellers-verified)", fontsize=9)
            axs[di][0].set_xlabel("edit distance d*"); axs[di][0].set_ylabel("segments")
        md.append("| k | segments with d* <= k | " + " | ".join(f"{a} ms" for a in ("Myers", "Edlib", "WFA2", "filter (pure)", "filter_fa (pure)")) +
                  " | filter hits | filter_fa hits | scans chosen by the fallback rule (filter) | one-hit read recall | two-hit read recall | one-hit end recall | two-hit end recall |")
        md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        rec1, rec2 = [], []
        for k in KS:
            res = {}
            truth = None
            for a in ALGOS:
                if a in ("edlib", "wfa"):
                    t0 = time.perf_counter()
                    for sg in segs:
                        if a == "edlib":
                            edlib.align(sg, t, mode="HW", task="locations", k=k)
                        else:
                            pywfa.WavefrontAligner(sg, span="ends-free", pattern_begin_free=0, pattern_end_free=0,
                                                   text_begin_free=len(t), text_end_free=len(t), distance="levenshtein")(t)
                    res[a] = dict(ms=(time.perf_counter() - t0) * 1e3 / NSEG)
                    continue
                ba, nofb = ("filter", True) if a == "filter_nofb" else (("filter_fa", True) if a == "filter_fa_nofb" else (a, False))
                info, sets = run_bin(ba, tpath, segs, k, nofallback=nofb, tag="dv")
                if a == "sellers":
                    truth = sets
                res[a] = dict(ms=info["time_ms"] / NSEG, hits=info["hits"] / NSEG, fb=int(info["fallbacks"]), sets=sets)
            mappable = [i for i in range(NSEG) if d[i] <= k]
            for a in EXACT:
                ok = all(res[a]["sets"].get(i, set()) == truth.get(i, set()) for i in range(NSEG))
                if not ok:
                    raise SystemExit(f"EXACT MISMATCH in divergence experiment: {a} k={k}")
            tot = sum(len(truth.get(i, ())) for i in mappable) or 1
            rr, re_ = {}, {}
            for a in ("blast1", "blast2"):
                st = res[a]["sets"]
                assert all(st.get(i, set()) <= truth.get(i, set()) for i in range(NSEG)), a + " unsound"
                rr[a] = sum(1 for i in mappable if st.get(i)) / max(1, len(mappable))
                re_[a] = sum(len(st.get(i, set()) & truth.get(i, set())) for i in mappable) / tot
            rec1.append(rr["blast1"]); rec2.append(rr["blast2"])
            rows.append(dict(direction=f"{qname}->{tname}", m=m, k=k, mappable=len(mappable), **{f"{a}_ms": res[a]["ms"] for a in ALGOS},
                             filter_hits=res["filter_nofb"]["hits"], filter_fa_hits=res["filter_fa_nofb"]["hits"],
                             recall_one_read=rr["blast1"], recall_two_read=rr["blast2"], recall_one_end=re_["blast1"], recall_two_end=re_["blast2"]))
            md.append(f"| {k} | {len(mappable)} | {res['myers']['ms']:.4f} | {res['edlib']['ms']:.4f} | {res['wfa']['ms']:.4f} | "
                      f"{res['filter_nofb']['ms']:.4f} | {res['filter_fa_nofb']['ms']:.4f} | {res['filter_nofb']['hits']:.0f} | {res['filter_fa_nofb']['hits']:.0f} | "
                      f"{res['filter']['fb']} / {NSEG} | {rr['blast1']:.2f} | {rr['blast2']:.2f} | {re_['blast1']:.2f} | {re_['blast2']:.2f} |")
        if m == 100:
            axs[di][1].plot(KS, rec1, "x-", color="tab:red", label="one-hit read recall"); axs[di][1].plot(KS, rec2, "+-", color="salmon", label="two-hit read recall")
            axs[di][1].set_ylim(0, 1.05); axs[di][1].set_xlabel("k"); axs[di][1].set_ylabel("recall of mappable segments (d* <= k)")
            axs[di][1].set_title(f"Heuristic recall, {qname} vs {tname}, real homologs", fontsize=9); axs[di][1].grid(alpha=.3); axs[di][1].legend(fontsize=8)
md.append("\nTimes are mean ms per segment on a 16.5 kbp text (all exact methods agree with Sellers on every segment and k; Edlib and WFA2 "
          "return the best match per read and include Python overhead). With a text this small the index costs about 20 ms to build "
          "(clearing the 4^11 table), so the filter is only worthwhile over many reads. \"pure\" = fallback disabled; the fallback column counts how many segments the cost-model rule sent to a full Myers scan.")
fig.tight_layout(); fig.savefig(os.path.join(RES, "divergence.png"), dpi=160); plt.close(fig)
with open(os.path.join(RES, "divergence.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
open(os.path.join(RES, "divergence.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md[:14]))
