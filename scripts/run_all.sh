#!/usr/bin/env bash
# One command: build, fetch data (skipped if present), crosscheck, calibrate cost model (NNLS), instruction profile, benchmark
# (headline experiment with 100 reads per setting), stress, plots, summary.
#   scripts/run_all.sh            full run (about 15 minutes)
#   scripts/run_all.sh --quick    smoke run (about 2 minutes, tiny sample sizes)
set -euo pipefail
cd "$(dirname "$0")/.."
QUICK="${1:-}"
echo "== build"
g++ -O3 -march=native -std=c++17 -Wall -o approx src/approx.cpp
python3 -c "import numpy, scipy, matplotlib, edlib, pywfa" 2>/dev/null || pip install -q numpy scipy matplotlib edlib pywfa
echo "== data"
if [ ! -s data/genome_synth.txt ]; then python3 scripts/simulate.py --n 4600000 --out data/genome_synth.txt; fi
if [ ! -s data/ecoli536.txt ] || [ ! -s data/genome.txt ]; then python3 scripts/fetch_data.py || echo "fetch failed: synthetic fallback"; fi
if [ ! -s data/genome.txt ]; then cp data/genome_synth.txt data/genome.txt; fi
head -3 data/SOURCES.txt
echo "== correctness (all exact algorithms vs Sellers; Edlib vs min-distance ends)"
python3 tests/edlib_semantics.py
python3 tests/crosscheck.py
APPROX_NOFALLBACK=1 python3 tests/crosscheck.py     # pure filter path, no scan fallback
echo "== cost-model calibration (synthetic genome; evaluated later on E. coli 536)"
python3 scripts/calibrate.py > results/calibration_log.txt
echo "== instructions per cell (cachegrind; perf hardware counters unavailable in this VM)"
python3 scripts/instr_profile.py > results/instr_profile_log.txt || echo "instr_profile failed"
echo "== benchmarks + stress"
python3 scripts/run_bench.py $QUICK | tee results/bench_log.txt | tail -3
echo "== plots + summary"
python3 scripts/plot.py
echo "done: results/summary.md, results/*.png, results/bench.csv"
