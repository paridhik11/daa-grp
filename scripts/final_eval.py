#!/usr/bin/env python3
"""Final pass: interleaved per-read timing, 3 repeats, ratios to Myers.  One process per setting; all methods run on every read in a rotating
order (`approx multi:...`), index build and startup excluded.  Speedup of a method on a read = time(Myers scan) / time(method), same read, same repeat;
tables give the median and the 10th to 90th percentile over reads x repeats.  Parameters are frozen (nothing here changes the algorithms).
E. coli 536 is the text the constants were tuned on (not blind); B. fragilis was never used for tuning (BLIND).
Writes results/final_eval.md/.csv.  Usage: [--quick]"""
import csv, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, BIN  # noqa: E402
from simulate import make_reads  # noqa: E402

QUICK = "--quick" in sys.argv
NR, REPS = (20, 1) if QUICK else (60, 3)
METHODS = ["myers", "lv", "filter_fa", "filter_fa_wg", "filter_s2", "filter_s3", "filter_fs2", "filter_fs3", "filter_auto"]
SET = [(100, k) for k in ([2, 8, 12, 17, 20] if QUICK else [1, 2, 4, 6, 8, 10, 12, 15, 17, 20, 25])] + ([] if QUICK else [(200, 30), (200, 40)])
rows = []


def one(text_path, label, m, k, seed):
    t = open(text_path).read().strip()
    reads = [r[0] for r in make_reads(t, NR, m, k, seed=seed)]
    pf = os.path.join(DATA, "_fe_reads.txt")
    open(pf, "w").write("\n".join(reads) + "\n")
    env = dict(os.environ, APPROX_REPS=str(REPS), APPROX_NOFALLBACK="0")
    r = subprocess.run([BIN, "multi:" + ",".join(METHODS), text_path, pf, str(k)], capture_output=True, text=True, check=True, env=env)
    ns = {}
    mism = None
    for l in r.stdout.splitlines():
        p = l.split()
        if p[0] == "I":
            ns[(p[1], int(p[2]), int(p[3]))] = float(p[4])
        elif p[0] == "MULTI":
            mism = int(p[4].split("=")[1])
    assert mism == 0, f"exactness mismatch {label} m={m} k={k}"
    out = {}
    for meth in METHODS[1:]:
        sp = np.array([ns[("myers", rp, i)] / ns[(meth, rp, i)] for rp in range(REPS) for i in range(NR)])
        out[meth] = (np.median(sp), np.percentile(sp, 10), np.percentile(sp, 90))
    my_ms = np.median([ns[("myers", rp, i)] for rp in range(REPS) for i in range(NR)]) / 1e6
    rows.append((label, m, k, my_ms, out))
    print(label, m, k, "ok", flush=True)


for label, f, seed in (("E. coli 536 (tuned on)", "ecoli536.txt", 101), ("B. fragilis (BLIND)", "bfragilis.txt", 202)):
    for m, k in SET:
        one(os.path.join(DATA, f), label, m, k, seed + k)

md = ["## Final pass: interleaved per-read timing, speedup over Myers\n",
      f"Every method runs on every read inside one process, in a rotating order; {REPS} repeats of {NR} reads per setting; index build and process startup are excluded. "
      "Cell = median speedup over Myers (time(Myers)/time(method) on the same read), with the 10th to 90th percentile in brackets (>1 faster than the scan). All methods return exactly the same end positions in every run (checked per read). "
      "Fallback active. Parameters frozen. The blind table is the B. fragilis one.\n"]
for label in ("E. coli 536 (tuned on)", "B. fragilis (BLIND)"):
    md += [f"\n### {label}\n", "| m, k | " + " | ".join(METHODS[1:]) + " |", "|---|" + "---|" * (len(METHODS) - 1)]
    for lb, m, k, my, out in rows:
        if lb != label:
            continue
        md.append(f"| m={m}, k={k} | " + " | ".join(f"{out[x][0]:.2f} ({out[x][1]:.2f} to {out[x][2]:.2f})" for x in METHODS[1:]) + " |")
    md.append("\nMedian absolute time of the Myers scan (ms per read, for scale only): " + ", ".join(f"m={m},k={k}: {my:.1f}" for lb, m, k, my, _ in rows if lb == label) + ".\n")
open(os.path.join(RES, "final_eval.md"), "w").write("\n".join(md) + "\n")
with open(os.path.join(RES, "final_eval.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["text", "m", "k", "myers_ms"] + [f"{x}_{q}" for x in METHODS[1:] for q in ("med", "p10", "p90")])
    for lb, m, k, my, out in rows:
        w.writerow([lb, m, k, f"{my:.3f}"] + [f"{out[x][i]:.3f}" for x in METHODS[1:] for i in range(3)])
