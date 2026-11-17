#!/usr/bin/env python3
"""Peak memory per method (measured): peak resident set size (VmHWM) of the process on E. coli 536, m = 100, k = 10, 5 reads.
Over baseline = VmHWM minus the steady VmRSS after loading the text (the memory the method itself needs).  Edlib and WFA2 are measured in a
separate Python process with ru_maxrss relative to a Python process that only loaded the text.  Writes results/memory.md and memory.csv."""
import csv
import os
import resource
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin  # noqa: E402
from simulate import make_reads  # noqa: E402

tp = os.path.join(DATA, "ecoli536.txt")
g = open(tp).read().strip()
n = len(g)
k = 10
reads = [r[0] for r in make_reads(g, 5, 100, k, seed=7)]
theory = {"noop": "text only (n bytes) plus read buffers", "sellers": "O(m)", "banded": "O(m)", "myers": "O(B sigma)",
          "myers_pruned": "O(B sigma)", "lv": "packed text n/4 bytes (1.2 MB) + L1-sized chunk buffers", "lv_direct": "2 layers of n + 2k ints plus n flags (about 44 MB)", "filter_auto": "same index",
          "filter": "index 4(4^11+1) + 4n bytes (37 MB) plus build temporaries (counter copy 17 MB, flags 5 MB)", "filter_fa": "same index", "filter_fa_w": "same index",
          "blast1": "same index", "blast2": "same index", "map": "same index + O(m (m + k)) per cluster"}
rows = []
base = None
for a in ("noop", "sellers", "banded", "myers", "myers_pruned", "lv", "lv_direct", "filter", "filter_auto", "filter_fa", "filter_fa_w", "blast1", "blast2", "map"):
    info, _ = run_bin(a, tp, reads, k, dump=False, tag="mem")
    kb = info["peak_kb"]
    if a == "noop":
        base = info["loaded_kb"]; noop_peak = kb
    over = (kb - info["loaded_kb"]) / 1024.0 if kb > noop_peak + 1024 else None   # None: hidden by the transient loading peak
    rows.append(dict(algo=a, peak_mb=kb / 1024.0, over_baseline_mb=over, theory=theory[a]))
code = r'''
import resource, sys
sys.path.insert(0, "%s")
txt = open("%s").read().strip().encode()
base = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
rd = open("%s").read().split()
if sys.argv[1] == "edlib":
    import edlib
    for r in rd: edlib.align(r, txt, mode="HW", task="locations", k=%d)
else:
    import pywfa
    for r in rd: pywfa.WavefrontAligner(r, span="ends-free", pattern_begin_free=0, pattern_end_free=0, text_begin_free=len(txt), text_end_free=len(txt), distance="levenshtein")(txt)
print(base, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
''' % (os.path.dirname(os.path.abspath(__file__)), tp, os.path.join(DATA, "_mem_reads.txt"), k)
open(os.path.join(DATA, "_mem_reads.txt"), "w").write("\n".join(reads) + "\n")
for a in ("edlib", "wfa"):
    out = subprocess.run([sys.executable, "-c", code, a], capture_output=True, text=True).stdout.split()
    if len(out) == 2:
        b0, b1 = int(out[0]), int(out[1])
        rows.append(dict(algo=a + " (Python process)", peak_mb=b1 / 1024.0, over_baseline_mb=(b1 - b0) / 1024.0 if b1 > b0 + 1024 else None,
                         theory="library; measured as ru_maxrss increase over loading the text"))
md = ["## Peak memory per method (measured)\n",
      f"E. coli 536 (n = {n:,}), m = 100, k = 10, 5 reads; peak resident set size (VmHWM). working memory = peak RSS minus the steady RSS after the text and reads are loaded and the heap trimmed ({base / 1024:.1f} MB); methods whose peak does not exceed the transient peak of loading the file (about 20 MB, the file is read line by line into a string) are shown as negligible. "
      "Edlib and WFA2 run in a Python process, so their baseline differs.\n",
      "| method | peak RSS (MB) | working memory of the method (MB) | expected |", "|---|---|---|---|"]
for r in rows:
    ob = "negligible (below the loading peak)" if r["over_baseline_mb"] is None else f"{r['over_baseline_mb']:.1f}"
    md.append(f"| {r['algo']} | {r['peak_mb']:.1f} | {ob} | {r['theory']} |")
open(os.path.join(RES, "memory.md"), "w").write("\n".join(md) + "\n")
with open(os.path.join(RES, "memory.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("\n".join(md))
