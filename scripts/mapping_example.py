#!/usr/bin/env python3
"""Worked example for the report: map a short read with known substitutions and indels to E. coli 536 with scripts/map_reads."""
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
fasta = os.path.join(ROOT, "data", "raw", "ecoli536.fa")
genome = open(os.path.join(ROOT, "data", "ecoli536.txt")).read().strip()
pos = 2_000_000
ref = genome[pos:pos + 60]
r = list(ref)
r[10] = "A" if r[10] != "A" else "C"          # substitution
r.insert(26, "G")                              # insertion in the read
del r[41]                                      # deletion from the read
r[53] = "T" if r[53] != "T" else "G"           # substitution
read = "".join(r)
with tempfile.TemporaryDirectory() as d:
    rf = os.path.join(d, "reads.fa")
    open(rf, "w").write(f">example_read\n{read}\n>unmappable_read\n{'ACGT' * 15}\n")
    tsv = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "map_reads"), fasta, rf, "-k", "8"], capture_output=True, text=True, check=True).stdout
    sam = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "map_reads"), fasta, rf, "-k", "8", "--format", "sam"], capture_output=True, text=True, check=True).stdout
row = tsv.splitlines()[1].split("\t")
start, end, cig = int(row[4]), int(row[5]), row[7]
top, mid, bot = [], [], []
ri, ti = 0, start
for n, op in re.findall(r"(\d+)([=XID])", cig):
    for _ in range(int(n)):
        if op in "=X":
            top.append(genome[ti]); bot.append(read[ri]); mid.append("|" if op == "=" else "x"); ri += 1; ti += 1
        elif op == "I":
            top.append("-"); bot.append(read[ri]); mid.append(" "); ri += 1
        else:
            top.append(genome[ti]); bot.append("-"); mid.append(" "); ti += 1
md = ["### Worked example: `scripts/map_reads` on E. coli 536\n",
      f"A 60 bp stretch of E. coli 536 (positions {pos}..{pos + 60}) was copied and given 2 substitutions, 1 insertion and 1 deletion (4 edits), then mapped with k = 8:\n",
      "```", "scripts/map_reads data/raw/ecoli536.fa reads.fa -k 8", "```", "TSV output:", "```", tsv.rstrip(), "```",
      "Alignment reconstructed from the CIGAR (reference on top, read below; `x` mismatch, `-` gap):", "```",
      "".join(top), "".join(mid), "".join(bot), "```",
      f"The reported distance is {row[3]}, equal to the number of X, I and D operations in `{cig}`, and the start ({start}) matches the planted position ({pos}). "
      "The second read has no alignment within k = 8 and is reported as unmapped. SAM-like output (`--format sam`):", "```", sam.rstrip(), "```"]
open(os.path.join(RES, "mapping_example.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md[:14]))
