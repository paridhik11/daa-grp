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
echo "== real-data cross-check of every exact algorithm vs Sellers"
python3 tests/crosscheck_real.py
python3 tests/test_pruned.py
python3 tests/test_lv.py
python3 tests/test_filter_s.py
echo "== alignment output validation (CIGAR) + worked example"
python3 tests/test_mapping.py
python3 scripts/mapping_example.py
echo "== cost-model calibration (synthetic genome; evaluated later on E. coli 536)"
python3 scripts/calibrate.py > results/calibration_log.txt
echo "== instructions per cell (cachegrind; perf hardware counters unavailable in this VM)"
python3 scripts/instr_profile.py > results/instr_profile_log.txt || echo "instr_profile failed"
echo "== benchmarks + stress"
python3 scripts/run_bench.py $QUICK | tee results/bench_log.txt | tail -3
echo "== seed-selection objective test (hits vs hits + lambda * windows)"
python3 scripts/objective_test.py > results/objective_test_log.txt
echo "== expected-case analysis (formula vs measured hits)"
python3 scripts/expected_hits.py > results/expected_hits_log.txt
echo "== real homologous sequences (mt_orang vs mt_human)"
python3 scripts/real_divergence.py > results/divergence_log.txt
echo "== phase 3A: filter_auto and cheaper selection, packed Landau-Vishkin, single-word Myers (targeted runs)"
python3 scripts/phase3a_filter.py > results/phase3a_filter_log.txt
python3 scripts/phase3a_noise.py >> results/phase3a_filter_log.txt
python3 scripts/phase3a_lv.py > results/phase3a_lv_log.txt
python3 scripts/phase3a_myers.py > results/phase3a_myers_log.txt
python3 scripts/phase3b_filter_s.py > results/phase3b_filter_s_log.txt
python3 scripts/phase3b_gate.py > results/phase3b_gate_log.txt
python3 scripts/phase3b_auto.py > results/phase3b_auto_log.txt
echo "== peak memory per method; DP-matrix figure"
python3 scripts/memory.py > results/memory_log.txt
python3 scripts/fig_dp_matrix.py
echo "== plots + summary"
python3 scripts/plot.py
echo "done: results/summary.md, results/*.png, results/bench.csv"
