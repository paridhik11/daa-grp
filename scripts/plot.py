#!/usr/bin/env python3
"""results/bench.csv -> results/*.png and results/summary.md.   (Every plot title names its dataset.)"""
import csv
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
rows = list(csv.DictReader(open(os.path.join(RES, "bench.csv"))))
CM = json.load(open(os.path.join(RES, "costmodel.json")))
NAME = {"ecoli536": "E. coli 536", "sarscov2": "SARS-CoV-2 (MT192765.1)", "lambda": "phage lambda",
        "mt_human": "human mtDNA", "genome_synth": "synthetic genome"}
MAIN = rows[0]["dataset"] if rows else "ecoli536"
MAINL = NAME.get(MAIN, MAIN)
LABEL = {"sellers": "Sellers DP", "lv": "Landau-Vishkin", "banded": "Ukkonen banded", "myers_pruned": "Myers + block pruning", "myers": "Myers bit-vector",
         "edlib": "Edlib (best match per read)", "wfa": "WFA2 (best match per read)",
         "filter": "filter, equal split", "filter_nofb": "filter, equal split (no fallback)",
         "filter_fa": "filter_fa, frequency-aware", "filter_fa_nofb": "filter_fa (no fallback)",
         "blast1": "BLAST-style one-hit", "blast2": "BLAST-style two-hit"}
STYLE = {"sellers": ("tab:gray", "o", "-"), "lv": ("tab:olive", "h", "-"), "banded": ("tab:orange", "s", "-"), "myers": ("tab:blue", "^", "-"), "myers_pruned": ("tab:cyan", "v", "-"),
         "edlib": ("tab:purple", "v", "--"), "wfa": ("tab:brown", "<", "--"),
         "filter": ("tab:green", "D", "-"), "filter_nofb": ("tab:green", "D", ":"),
         "filter_fa": ("darkgreen", "P", "-"), "filter_fa_nofb": ("darkgreen", "P", ":"),
         "blast1": ("tab:red", "x", "-"), "blast2": ("salmon", "+", "-")}


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def sel(exp, algo, **kw):
    r = [x for x in rows if x["exp"] == exp and x["algo"] == algo and all(x[k] == str(v) for k, v in kw.items())]
    try:
        return sorted(r, key=lambda x: f(x["x"]))
    except Exception:
        return r


def series(exp, algo, key, **kw):
    r = sel(exp, algo, **kw)
    return [f(x["x"]) for x in r], [f(x[key]) for x in r]


def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(RES, name), dpi=160); plt.close(fig)


def crossover(xs, a, b):
    """first x where a >= b, log-linear interpolation of the ratio; None if never."""
    pts = [(x, math.log(max(u, 1e-12) / max(v, 1e-12))) for x, u, v in zip(xs, a, b)]
    for (x0, r0), (x1, r1) in zip(pts, pts[1:]):
        if r0 < 0 <= r1:
            return x0 + (x1 - x0) * (0 - r0) / (r1 - r0)
    return None if not pts or pts[-1][1] < 0 else pts[0][0]


md = ["# Benchmark summary\n",
      f"Main text: **{MAINL}** (see `data/SOURCES.txt`; no K-12 genome was reachable). Search time is per pattern, "
      "single thread; index build is reported separately.\n",
      "Edlib and WFA2 return the BEST match per read (minimum distance, ends at that minimum), not all matches "
      "within k, and include Python call overhead; compare them on time only.\n",
      f"Cost-model constants (ns), calibrated on {CM.get('calibrated_on')} and evaluated on {CM.get('evaluated_on')}: "
      f"c_hit={CM['c_hit']:.1f}, c_win={CM['c_win']:.0f}, c_bp={CM['c_bp']:.2f}, a_scan={CM['a_scan']:.2f}, "
      f"c_sel={CM['c_sel']:.1f}\n"]


# ---------------------------------------------------------------- headline experiment (100 reads): median and 10th-90th percentile
HL = {}
for r in csv.DictReader(open(os.path.join(RES, "headline.csv"))):
    HL[(r["algo"], int(r["k"]))] = {k_: (v if k_ == "algo" else float(v)) for k_, v in r.items()}
HK = sorted({k_ for _, k_ in HL})
NR = int(HL[("myers", HK[0])]["npat"])


def hv(a, k, key):
    return HL[(a, k)][key]


def band(a, key):  # (xs, med, yerr) with asymmetric 10th-90th percentile error bars
    pre = {"t": "t", "hits": "hits", "bp": "bp"}[key]
    xs_ = [k_ / 100 for k_ in HK]
    med = [hv(a, k_, pre + "_med") for k_ in HK]
    lo = [max(0.0, m_ - hv(a, k_, pre + "_p10")) for k_, m_ in zip(HK, med)]
    hi = [max(0.0, hv(a, k_, pre + "_p90") - m_) for k_, m_ in zip(HK, med)]
    return xs_, med, [lo, hi]


def mp(a, k, key, fmt):  # "median [p10, p90]"
    pre = {"t": "t", "hits": "hits", "bp": "bp"}[key]
    return f"{fmt % hv(a, k, pre + '_med')} [{fmt % hv(a, k, pre + '_p10')}, {fmt % hv(a, k, pre + '_p90')}]"


md.append(f"\nHeadline experiment ({MAINL}, m = 100, {NR} reads per setting): values are the **median [10th, 90th percentile] over reads** "
          f"of per-read query time; other experiments use means over 20 reads (Sellers 1, WFA2 5).\n")

# ---------------------------------------------------------------- 1. time vs error (all methods)
fig, ax = plt.subplots(figsize=(8, 5))
for a in ["sellers", "lv", "banded", "myers_pruned", "edlib", "wfa", "blast1", "blast2"]:
    xs, ys = series("error", a, "ms_per_pattern")
    c_, mk, ls = STYLE[a]
    ax.plot(xs, ys, color=c_, marker=mk, ls=ls, label=LABEL[a], ms=5)
for a in ["myers", "filter", "filter_fa"]:
    xs, med, yerr = band(a, "t")
    c_, mk, ls = STYLE[a]
    ax.errorbar(xs, med, yerr=yerr, color=c_, marker=mk, ls=ls, label=LABEL[a] + " (median, 10-90%)", ms=5, capsize=2, lw=1.4)
ax.set_yscale("log"); ax.set_xlabel("error rate k/m (m = 100)"); ax.set_ylabel("search time per pattern (ms, log)")
ax.set_title(f"Time vs error rate, {MAINL}"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7, ncol=2)
save(fig, "time_vs_error.png")

# ---------------------------------------------------------------- 2. filter_fa vs equal split (no fallback)
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
for ax, key, yl in zip(axs, ["hits", "bp", "t"],
                       ["seed hits (candidates) per pattern", "bp verified per pattern", "query time per pattern (ms)"]):
    for a in ("filter_nofb", "filter_fa_nofb"):
        xs, med, yerr = band(a, key)
        ax.errorbar(xs, med, yerr=yerr, color=STYLE[a][0], marker=STYLE[a][1], label=LABEL[a], capsize=2, ms=5)
    ax.set_yscale("log"); ax.set_xlabel("error rate k/m (m = 100)"); ax.set_ylabel(yl + " (median, 10-90%)"); ax.grid(alpha=.3, which="both")
axs[0].legend(fontsize=8); fig.suptitle(f"Frequency-aware vs equal split (pure filter, no fallback), {MAINL}, {NR} reads")
save(fig, "fa_vs_equal.png")
md.append("\n## Frequency-aware vs equal split (pure filter, no fallback), " + MAINL + "\n")
md.append("| k | hits equal | hits fa | hit ratio (medians) | bp verified equal | bp verified fa | time equal (ms) | time fa (ms) | fa speedup vs equal | fa speedup vs Myers |")
md.append("|---|---|---|---|---|---|---|---|---|---|")
for k in HK:
    te, tf, tm = hv("filter_nofb", k, "t_med"), hv("filter_fa_nofb", k, "t_med"), hv("myers", k, "t_med")
    md.append(f"| {k} | {mp('filter_nofb', k, 'hits', '%.0f')} | {mp('filter_fa_nofb', k, 'hits', '%.0f')} | "
              f"{hv('filter_nofb', k, 'hits_med') / max(1, hv('filter_fa_nofb', k, 'hits_med')):.1f}x | "
              f"{mp('filter_nofb', k, 'bp', '%.0f')} | {mp('filter_fa_nofb', k, 'bp', '%.0f')} | "
              f"{mp('filter_nofb', k, 't', '%.3f')} | {mp('filter_fa_nofb', k, 't', '%.3f')} | {te / tf:.2f}x | {tm / tf:.0f}x |")

# ---------------------------------------------------------------- 3. cost model: predicted vs measured
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
xs = [k_ / 100 for k_ in HK]
meas_f = [hv("filter_nofb", k_, "t_med") for k_ in HK]; pred_f = [hv("filter_nofb", k_, "predf_med") for k_ in HK]
meas_fa = [hv("filter_fa_nofb", k_, "t_med") for k_ in HK]; pred_fa = [hv("filter_fa_nofb", k_, "predf_med") for k_ in HK]
meas_m = [hv("myers", k_, "t_med") for k_ in HK]; pred_s = [hv("filter_nofb", k_, "preds_med") for k_ in HK]
ax = axs[0]
for a, mk_, col, lab in (("filter_nofb", "D", "tab:green", "filter"), ("filter_fa_nofb", "P", "darkgreen", "filter_fa"), ("myers", "^", "tab:blue", "Myers")):
    xs_, med, yerr = band(a, "t")
    ax.errorbar(xs_, med, yerr=yerr, color=col, marker=mk_, ls="-", capsize=2, label=lab + " measured (median, 10-90%)")
ax.plot(xs, pred_f, "D:", color="tab:green", label="filter predicted"); ax.plot(xs, pred_fa, "P:", color="darkgreen", label="filter_fa predicted")
ax.plot(xs, pred_s, "^:", color="tab:blue", label="Myers predicted")
cm_m = crossover(xs, meas_f, meas_m); cm_p = crossover(xs, pred_f, pred_s)
for c_, col, ls, lab in ((cm_m, "k", "-", "measured crossover"), (cm_p, "gray", "--", "predicted crossover")):
    if c_ is not None:
        ax.axvline(c_, color=col, ls=ls, lw=1, label=f"{lab} k/m={c_:.3f}")
ax.set_yscale("log"); ax.set_xlabel("error rate k/m (m = 100)"); ax.set_ylabel("time per pattern (ms)")
ax.set_title(f"Cost model vs measurement, {MAINL}"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7)
ax = axs[1]
pm = [pf_ / m_ for pf_, m_ in zip(pred_f, meas_f)]; pma = [a_ / b_ for a_, b_ in zip(pred_fa, meas_fa)]; pms = [a_ / b_ for a_, b_ in zip(pred_s, meas_m)]
ax.plot(xs, pm, "o-", color="tab:green", label="filter predicted / measured")
ax.plot(xs, pma, "o-", color="darkgreen", label="filter_fa predicted / measured")
ax.plot(xs, pms, "o-", color="tab:blue", label="Myers predicted / measured")
ax.axhline(1, color="k", lw=.8); ax.set_yscale("log"); ax.set_xlabel("error rate k/m"); ax.set_ylabel("prediction ratio (median pred / median measured)")
ax.set_title("Prediction error (1 = exact)"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
save(fig, "costmodel.png")
cm_fa_m = crossover(xs, meas_fa, meas_m); cm_fa_p = crossover(xs, pred_fa, pred_s)
md.append("\n## Cost model (predicted vs measured), " + MAINL + "\n")
md.append(f"Fit: {CM.get('fit', 'least squares')} on {CM.get('calibrated_on')}; evaluated here on {MAINL} (a different text, but NOT blind: the calibration range was widened after a first fit failed on it; see the blind-evaluation section). "
          f"c_win = {CM['c_win']:.0f} ns (unconstrained least squares gave {CM.get('c_win_unconstrained', float('nan')):.0f}, so the sign constraint "
          f"{'was not' if CM['c_win'] > 0 else 'was'} active). Hit counts come from the index; windows and merged bp are estimated assuming uniform hit placement.\n")
fmt = lambda c_: "never within tested k" if c_ is None else f"k/m = {c_:.3f} (k = {c_ * 100:.1f})"
md.append(f"- Crossover filter vs Myers: measured {fmt(cm_m)}, predicted {fmt(cm_p)}")
md.append(f"- Crossover filter_fa vs Myers: measured {fmt(cm_fa_m)}, predicted {fmt(cm_fa_p)}")
for name, ratios, meas in (("filter", pm, meas_f), ("filter_fa", pma, meas_fa)):
    good = [r_ for r_, mv in zip(ratios, meas) if mv > 0.01]
    lr = sorted(abs(math.log(r_)) for r_ in good)
    md.append(f"- {name} time prediction: median factor error {math.exp(lr[len(lr) // 2]):.2f}x, worst {math.exp(lr[-1]):.2f}x over {len(good)} values of k "
              f"(measured median > 0.01 ms; prediction ratio uses medians)")
lrm = sorted(abs(math.log(r_)) for r_ in pms)
md.append(f"- Myers scan prediction: median factor error {math.exp(lrm[len(lrm) // 2]):.2f}x, worst {math.exp(lrm[-1]):.2f}x")
md.append("\n| k | filter measured (ms) | predicted | filter_fa measured | predicted | Myers measured | predicted |")
md.append("|---|---|---|---|---|---|---|")
for i, k in enumerate(HK):
    md.append(f"| {k} | {mp('filter_nofb', k, 't', '%.3f')} | {pred_f[i]:.3f} | {mp('filter_fa_nofb', k, 't', '%.3f')} | {pred_fa[i]:.3f} | "
              f"{mp('myers', k, 't', '%.1f')} | {pred_s[i]:.1f} |")

# ---------------------------------------------------------------- 4. break-even
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), sharey=False)
md.append("\n## Index build vs query time and break-even, " + MAINL + "\n")
md.append("Q* = index_ms / (Myers median ms per query - filter median ms per query); the number of queries after which building "
          "the index plus filtering beats repeated plain Myers scans (never, if the filter is not faster). Query times are medians over "
          f"{NR} reads (with the cost-model fallback active); the index build time is one measurement per setting.\n")
md.append("| k | index build (ms) | Myers (ms/query) | filter (ms/query) | Q* filter | filter_fa (ms/query) | Q* filter_fa |")
md.append("|---|---|---|---|---|---|---|")
def qstar(idx, tm, tf):
    return idx / (tm - tf) if tm > tf else float("inf")
for k in HK:
    tm, t1, t2 = hv("myers", k, "t_med"), hv("filter", k, "t_med"), hv("filter_fa", k, "t_med")
    md.append(f"| {k} | {hv('filter', k, 'index_ms'):.0f} | {tm:.1f} | {t1:.3f} | {qstar(hv('filter', k, 'index_ms'), tm, t1):.2f} | "
              f"{t2:.3f} | {qstar(hv('filter_fa', k, 'index_ms'), tm, t2):.2f} |")
for ax, k in zip(axs, [k_ for k_ in (4, 8, 12) if k_ in HK]):
    Q = list(range(0, 21))
    tm, t1, t2 = hv("myers", k, "t_med"), hv("filter", k, "t_med"), hv("filter_fa", k, "t_med")
    i1, i2 = hv("filter", k, "index_ms"), hv("filter_fa", k, "index_ms")
    ax.plot(Q, [q * tm for q in Q], "^-", color="tab:blue", label="Myers (no index)")
    ax.plot(Q, [i1 + q * t1 for q in Q], "D-", color="tab:green", label="filter + index")
    ax.plot(Q, [i2 + q * t2 for q in Q], "P-", color="darkgreen", label="filter_fa + index")
    if tm > t1:
        qs = qstar(i1, tm, t1)
        ax.axvline(qs, color="k", ls="--", lw=1); ax.text(qs + .3, i1 * .5, f"Q* = {qs:.1f}", fontsize=8)
    ax.set_title(f"k = {k} (m = 100), {MAINL}"); ax.set_xlabel("number of queries"); ax.set_ylabel("cumulative ms")
    ax.grid(alpha=.3); ax.legend(fontsize=7)
save(fig, "breakeven.png")

# ---------------------------------------------------------------- 4b. blind evaluation of the cost model on three texts
BLIND = [("synth", "Synthetic genome X (fitted text, new reads; in-sample)", "headline_synth_insample.csv"),
         ("ecoli", "E. coli 536 (not fitted; calibration range chosen after seeing it; NOT blind)", "headline.csv"),
         ("bfrag", "B. fragilis (never seen before the constants were frozen; BLIND)", "headline_bfragilis.csv")]
blind_rows = {}
for key_, lab_, fn_ in BLIND:
    path_ = os.path.join(RES, fn_)
    if not os.path.exists(path_):
        continue
    d_ = {}
    for r in csv.DictReader(open(path_)):
        d_[(r["algo"], int(r["k"]))] = {k2: (v if k2 == "algo" else float(v)) for k2, v in r.items()}
    blind_rows[key_] = (lab_, d_, sorted({k2 for _, k2 in d_}))
if blind_rows:
    fig, axs = plt.subplots(1, len(blind_rows), figsize=(5.4 * len(blind_rows), 4.4), squeeze=False)
    md.append("\n## Cost model: which figures are blind and which are fitted\n")
    md.append(f"Protocol (fixed before the blind text was read): constants are fitted by NNLS on the synthetic genome X ({CM.get('calibrated_on')}) with k up to 25, "
              "then frozen. Prediction error = median over k of the factor max(pred/measured, measured/pred) between the median predicted and the median measured "
              "per-read time of the pure filter (fallback off), for k where the measured median exceeds 0.01 ms. The fit was not changed after seeing any evaluation text.\n")
    md.append("| text | status | reads per k | filter error (median / worst) | filter_fa error (median / worst) | Myers error (median) | crossover filter measured / predicted (k/m) | crossover filter_fa measured / predicted (k/m) |")
    md.append("|---|---|---|---|---|---|---|---|")
    for ax, (key_, (lab_, d_, ks_)) in zip(axs[0], blind_rows.items()):
        g_ = lambda a, k_, f_: d_[(a, k_)][f_]
        xs_ = [k_ / 100 for k_ in ks_]
        mf = [g_("filter_nofb", k_, "t_med") for k_ in ks_]; pf_ = [g_("filter_nofb", k_, "predf_med") for k_ in ks_]
        mfa = [g_("filter_fa_nofb", k_, "t_med") for k_ in ks_]; pfa = [g_("filter_fa_nofb", k_, "predf_med") for k_ in ks_]
        mm = [g_("myers", k_, "t_med") for k_ in ks_]; ps_ = [g_("filter_nofb", k_, "preds_med") for k_ in ks_]
        ax.plot(xs_, mf, "D-", color="tab:green", label="filter measured"); ax.plot(xs_, pf_, "D:", color="tab:green", label="filter predicted")
        ax.plot(xs_, mfa, "P-", color="darkgreen", label="filter_fa measured"); ax.plot(xs_, pfa, "P:", color="darkgreen", label="filter_fa predicted")
        ax.plot(xs_, mm, "^-", color="tab:blue", label="Myers measured"); ax.plot(xs_, ps_, "^:", color="tab:blue", label="Myers predicted")
        cmx, cpx = crossover(xs_, mf, mm), crossover(xs_, pf_, ps_)
        cfx, cfp = crossover(xs_, mfa, mm), crossover(xs_, pfa, ps_)
        for c_, col_, ls_ in ((cmx, "k", "-"), (cpx, "gray", "--")):
            if c_ is not None:
                ax.axvline(c_, color=col_, ls=ls_, lw=1)
        ax.set_yscale("log"); ax.set_xlabel("error rate k/m (m = 100)"); ax.set_ylabel("time per pattern (ms)"); ax.grid(alpha=.3, which="both")
        ax.set_title(lab_, fontsize=8); ax.legend(fontsize=6)
        def err(meas, pred):
            r_ = sorted(max(p_ / m_, m_ / p_) for p_, m_ in zip(pred, meas) if m_ > 0.01)
            return (r_[len(r_) // 2], r_[-1]) if r_ else (float("nan"), float("nan"))
        ef, efa = err(mf, pf_), err(mfa, pfa)
        em = err(mm, ps_)
        fx = lambda c_: "never" if c_ is None else f"{c_:.3f}"
        status = {"synth": "fitted text, new reads", "ecoli": "not fitted, range chosen after seeing it", "bfrag": "blind"}[key_]
        md.append(f"| {lab_.split(' (')[0]} | {status} | {int(g_('myers', ks_[0], 'npat'))} | {ef[0]:.2f}x / {ef[1]:.2f}x | {efa[0]:.2f}x / {efa[1]:.2f}x | {em[0]:.2f}x | "
                  f"{fx(cmx)} / {fx(cpx)} | {fx(cfx)} / {fx(cfp)} |")
    save(fig, "costmodel_blind.png")
    md.append("\nFitted: the constants themselves and the synthetic-text row. Not blind: the E. coli row (the calibration range was widened to k <= 25 after "
              "a first fit restricted to k <= 15 failed on this text). Blind: the B. fragilis row only.\n")

# ---------------------------------------------------------------- 5. cells per second
fig, ax = plt.subplots(figsize=(7.5, 4.6))
md.append("\n## Cells updated per second, " + MAINL + " (m = 100)\n")
md.append("Sellers and Myers: m*n DP cells represented per pattern (Myers updates 64 per word operation); banded: cells "
          "actually computed.\n")
md.append("| k | Sellers cells/s | banded cells computed | banded cells/s | Myers cells represented | Myers cells/s | Myers/banded cells/s | banded ms | Myers ms |")
md.append("|---|---|---|---|---|---|---|---|---|")
B = {int(r["k"]): r for r in sel("error", "banded")}
M = {int(r["k"]): r for r in sel("error", "myers")}
S = {int(r["k"]): r for r in sel("error", "sellers")}
for a, lab in (("sellers", "Sellers DP"), ("banded", "Ukkonen banded"), ("myers", "Myers bit-vector")):
    xs_, ys_ = series("error", a, "cells_per_s")
    ax.plot(xs_, ys_, marker=STYLE[a][1], color=STYLE[a][0], label=lab)
for k in sorted(B):
    md.append(f"| {k} | {f(S[k]['cells_per_s']):.2e} | {f(B[k]['cells']):.2e} | {f(B[k]['cells_per_s']):.2e} | "
              f"{f(M[k]['cells']):.2e} | {f(M[k]['cells_per_s']):.2e} | {f(M[k]['cells_per_s']) / f(B[k]['cells_per_s']):.0f}x | "
              f"{f(B[k]['ms_per_pattern']):.1f} | {f(M[k]['ms_per_pattern']):.1f} |")
ax.set_yscale("log"); ax.set_xlabel("error rate k/m (m = 100)"); ax.set_ylabel("DP cells per second")
ax.set_title(f"Cell throughput, {MAINL}"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
save(fig, "cells_per_sec.png")

# ---------------------------------------------------------------- 6. stress
scen_order = ["uniform", "clustered", "indel_heavy", "lowcomplex_edge", "lowcomplex_inside"]
stress_k = sorted({int(r["k"]) for r in rows if r["exp"] == "stress"})
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
w = 0.35
for ai, key in enumerate(("recall_read", "recall_end")):
    ax = axs[ai]
    for ki, k in enumerate(stress_k):
        for hi, a in enumerate(("blast1", "blast2")):
            vals = []
            for sc in scen_order:
                r = [x for x in rows if x["exp"] == "stress" and x["algo"] == a and x["x"] == sc and int(x["k"]) == k]
                vals.append(f(r[0][key]) if r else float("nan"))
            pos = [i + (ki * 2 + hi - (len(stress_k) * 2 - 1) / 2) * (0.8 / (len(stress_k) * 2)) for i in range(len(scen_order))]
            ax.bar(pos, vals, width=0.8 / (len(stress_k) * 2), label=f"{LABEL[a]}, k={k}")
    ax.set_xticks(range(len(scen_order))); ax.set_xticklabels(scen_order, rotation=20, fontsize=8)
    ax.set_ylim(0, 1.05); ax.set_ylabel("recall of reads" if key == "recall_read" else "recall of match end positions")
    ax.set_title(f"Heuristic {'read' if ai == 0 else 'end-position'} recall, {MAINL} + low-complexity")
    ax.grid(alpha=.3, axis="y")
h_, l_ = axs[0].get_legend_handles_labels()
ax = axs[2]
for hi, a in enumerate(("filter_nofb", "filter_fa_nofb")):
    k = stress_k[min(1, len(stress_k) - 1)]
    vals = []
    for sc in scen_order:
        r = [x for x in rows if x["exp"] == "stress" and x["algo"] == a and x["x"] == sc and int(x["k"]) == k]
        vals.append(f(r[0]["hits_pp"]) if r else float("nan"))
    ax.bar([i + (hi - .5) * .38 for i in range(len(scen_order))], vals, width=.38, label=LABEL[a])
ax.set_yscale("log"); ax.set_xticks(range(len(scen_order))); ax.set_xticklabels(scen_order, rotation=20, fontsize=8)
ax.set_ylabel("seed hits per pattern"); ax.set_title(f"Filter candidates under stress (k={stress_k[min(1, len(stress_k) - 1)]})")
ax.legend(fontsize=7); ax.grid(alpha=.3, axis="y", which="both")
fig.legend(h_, l_, loc="lower center", ncol=6, fontsize=7)
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(os.path.join(RES, "stress.png"), dpi=160); plt.close(fig)
md.append("\n## Stress scenarios, " + MAINL + " + ~1% injected microsatellite/homopolymer arrays (m = 100, 20 reads each)\n")
md.append("Recall by read = fraction of reads (that have a true match within k) for which the heuristic reports at least one match; "
          "recall by end = fraction of true end positions reported. Exact methods are asserted equal to Myers in every row.\n")
md.append("| scenario | k | one-hit read | one-hit end | two-hit read | two-hit end | filter hits | filter_fa hits | filter ms | filter_fa ms | Myers ms |")
md.append("|---|---|---|---|---|---|---|---|---|---|---|")
for sc in scen_order:
    for k in stress_k:
        def g(a, key):
            r = [x for x in rows if x["exp"] == "stress" and x["algo"] == a and x["x"] == sc and int(x["k"]) == k]
            return r[0][key] if r else ""
        md.append(f"| {sc} | {k} | {f(g('blast1', 'recall_read')):.2f} | {f(g('blast1', 'recall_end')):.2f} | "
                  f"{f(g('blast2', 'recall_read')):.2f} | {f(g('blast2', 'recall_end')):.2f} | {f(g('filter_nofb', 'hits_pp')):.0f} | "
                  f"{f(g('filter_fa_nofb', 'hits_pp')):.0f} | {f(g('filter_nofb', 'ms_per_pattern')):.3f} | "
                  f"{f(g('filter_fa_nofb', 'ms_per_pattern')):.3f} | {f(g('myers', 'ms_per_pattern')):.1f} |")

# ---------------------------------------------------------------- 7. heuristic recall vs error, length, text size, datasets
fig, ax = plt.subplots(figsize=(6.6, 4.3))
for a in ("blast1", "blast2"):
    xs_, ys_ = series("error", a, "recall_read")
    ax.plot(xs_, ys_, marker=STYLE[a][1], color=STYLE[a][0], label=LABEL[a] + " (reads)")
    xs_, ys_ = series("error", a, "recall_end")
    ax.plot(xs_, ys_, marker=STYLE[a][1], color=STYLE[a][0], ls="--", label=LABEL[a] + " (end positions)")
ax.axhline(1.0, color="tab:green", label="all exact methods"); ax.set_ylim(0, 1.05)
ax.set_xlabel("error rate k/m (m = 100)"); ax.set_ylabel("recall"); ax.set_title(f"Heuristic recall vs error rate, {MAINL}")
ax.grid(alpha=.3); ax.legend(fontsize=7)
save(fig, "recall_vs_error.png")

for exp, xl, fn, ttl, logx in (("length", "pattern length m (5% error)", "time_vs_length.png", "Time vs pattern length", True),
                               ("textsize", "text length n (bp, prefixes)", "time_vs_n.png", "Time vs text size (m=100, k=5)", True)):
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for a in ["sellers", "lv", "banded", "myers", "myers_pruned", "edlib", "wfa", "filter", "filter_fa", "blast1", "blast2"]:
        xs_, ys_ = series(exp, a, "ms_per_pattern")
        if xs_:
            ax.plot(xs_, ys_, color=STYLE[a][0], marker=STYLE[a][1], ls=STYLE[a][2], label=LABEL[a], ms=5)
    ax.set_yscale("log"); ax.set_xlabel(xl); ax.set_ylabel("search time per pattern (ms, log)")
    if logx:
        ax.set_xscale("log")
    ax.set_title(f"{ttl}, {MAINL}"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7, ncol=2)
    save(fig, fn)

md.append("\n## Results by dataset (m = 100; ms per pattern; index build once)\n")
ds = [d for d in dict.fromkeys(r["dataset"] for r in rows if r["exp"] == "datasets")]
md.append("| dataset | length (bp) | k | Sellers | Myers | Edlib (best/read) | WFA2 (best/read) | filter | filter_fa | one-hit recall (read) | two-hit recall (read) | index build (ms) |")
md.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for d in ds:
    for k in sorted({int(r["k"]) for r in rows if r["exp"] == "datasets" and r["dataset"] == d}):
        def g(a, key):
            r = [x for x in rows if x["exp"] == "datasets" and x["dataset"] == d and x["algo"] == a and int(x["k"]) == k]
            return r[0][key] if r else ""
        md.append(f"| {NAME.get(d, d)} | {int(f(g('myers', 'n')))} | {k} | {f(g('sellers', 'ms_per_pattern')):.2f} | "
                  f"{f(g('myers', 'ms_per_pattern')):.2f} | {f(g('edlib', 'ms_per_pattern')):.2f} | {f(g('wfa', 'ms_per_pattern')):.2f} | "
                  f"{f(g('filter', 'ms_per_pattern')):.4f} | {f(g('filter_fa', 'ms_per_pattern')):.4f} | "
                  f"{f(g('blast1', 'recall_read')):.2f} | {f(g('blast2', 'recall_read')):.2f} | {f(g('filter', 'index_ms')):.1f} |")
fig, ax = plt.subplots(figsize=(8, 4.4))
labels = [f"{NAME.get(d, d)}\nk={k}" for d in ds for k in sorted({int(r['k']) for r in rows if r['exp'] == 'datasets' and r['dataset'] == d})]
combos = [(d, k) for d in ds for k in sorted({int(r['k']) for r in rows if r['exp'] == 'datasets' and r['dataset'] == d})]
for ai, a in enumerate(("myers", "edlib", "filter", "filter_fa")):
    vals = []
    for d, k in combos:
        r = [x for x in rows if x["exp"] == "datasets" and x["dataset"] == d and x["algo"] == a and int(x["k"]) == k]
        vals.append(f(r[0]["ms_per_pattern"]) if r else float("nan"))
    ax.bar([i + (ai - 1.5) * .2 for i in range(len(combos))], vals, width=.2, label=LABEL[a], color=STYLE[a][0])
ax.set_yscale("log"); ax.set_xticks(range(len(combos))); ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("ms per pattern (log)"); ax.set_title("Query time by dataset (m = 100)"); ax.legend(fontsize=7); ax.grid(alpha=.3, axis="y", which="both")
save(fig, "datasets.png")

# ---------------------------------------------------------------- headline numbers
md.append("\n## Key numbers (main text, m = 100)\n")
md.append(f"Myers, `filter`, `filter_fa`: median [10th, 90th percentile] over {NR} reads, ms per query. Other columns: mean over 20 reads (Sellers 1, WFA2 5); "
          "recall columns are heuristic read recall (20 reads).\n")
md.append("| k | k/m | Sellers | Banded | Myers | Edlib (best/read) | WFA2 (best/read) | filter | filter_fa | one-hit | two-hit | fallbacks (of 100) filter | filter_fa |")
md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for k in HK:
    def g(a, key):
        r_ = [x for x in rows if x["exp"] == "error" and x["algo"] == a and int(x["k"]) == k]
        return r_[0][key] if r_ else ""
    md.append(f"| {k} | {k / 100:.2f} | {f(g('sellers', 'ms_per_pattern')):.0f} | {f(g('banded', 'ms_per_pattern')):.1f} | "
              f"{mp('myers', k, 't', '%.1f')} | {f(g('edlib', 'ms_per_pattern')):.1f} | {f(g('wfa', 'ms_per_pattern')):.1f} | "
              f"{mp('filter', k, 't', '%.3f')} | {mp('filter_fa', k, 't', '%.3f')} | "
              f"{f(g('blast1', 'ms_per_pattern')):.3f} (rec {f(g('blast1', 'recall_read')):.2f}) | "
              f"{f(g('blast2', 'ms_per_pattern')):.3f} (rec {f(g('blast2', 'recall_read')):.2f}) | "
              f"{int(hv('filter', k, 'fallbacks'))} | {int(hv('filter_fa', k, 'fallbacks'))} |")
for extra_md in ("instr_profile.md", "expected_hits.md", "divergence.md", "mapping_validation.md", "mapping_example.md", "objective_test.md", "pruned_myers.md", "memory.md"):
    ipath = os.path.join(RES, extra_md)
    if os.path.exists(ipath):
        md.append("\n" + open(ipath).read())
open(os.path.join(RES, "summary.md"), "w").write("\n".join(md) + "\n")
print("plots and results/summary.md written")
