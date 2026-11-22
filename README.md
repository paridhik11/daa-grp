# Approximate DNA matching: comparative study (DAA project, Problem 8)

One command (build, fetch real genomes, correctness tests, NNLS cost-model calibration, instruction profile, benchmarks including a 100-read headline experiment with median and 10th-90th percentile, blind cost-model test on a second genome, stress tests, expected-case analysis, real-homolog experiment, memory, plots, summary; about 1 hour):

```
make bench
```

Smoke run: `make quick`. Correctness tests only: `make test`. Map reads with CIGAR output:

```
scripts/map_reads GENOME.fa READS.fa -k 8 [--format tsv|sam]
```

| Path | Contents |
|---|---|
| `report/REPORT.md` | Full report: problem, algorithms, proofs, complexity, conclusion, gaps, references |
| `results/summary.md` | All tables (break-even, cost model and blind test, cells per second, stress, datasets, expected-case analysis, real divergence, objective test, memory, mapping validation) |
| `results/*.png` | Plots (including `dp_matrix.png`, `costmodel_blind.png`, `expected_hits.png`); `results/bench.csv`, `headline*.csv` raw data |
| `src/approx.cpp` | Sellers, banded, Myers, Myers with block pruning, Landau-Vishkin (direct extension), equal-split filter, frequency-aware filter (hits, and hits + lambda * windows), `filter_auto` (per-query cost-model choice), BLAST-style heuristic, traceback (`map`) |
| `scripts/` | `fetch_data.py`, `simulate.py`, `calibrate.py`, `run_bench.py`, `plot.py`, `phase3b_*.py`, `final_eval.py`, `final_gate.py`, `instr_profile.py`, `expected_hits.py`, `real_divergence.py`, `objective_test.py`, `memory.py`, `fig_dp_matrix.py`, `mapping_example.py`, `phase3a_filter.py`, `phase3a_noise.py`, `phase3a_lv.py`, `phase3a_myers.py`, `map_reads`, `run_all.sh` |
| `tests/` | `crosscheck.py`, `crosscheck_real.py`, `test_pruned.py`, `test_lv.py`, `test_mapping.py` (all against Sellers), `edlib_semantics.py` |
| `data/SOURCES.txt` | Dataset provenance (URL, header, length, GC, SHA-256) |

Datasets: E. coli 536 (main), B. fragilis (blind evaluation only), SARS-CoV-2, phage lambda, human and orangutan mtDNA, fetched from GitHub raw files. K-12 MG1655 could not be obtained (NCBI, ENA, UCSC and Ensembl were unreachable), so no dataset is called K-12. Edlib and WFA2 return the best match per read, not all matches within k. Hardware performance counters are unavailable in the sandbox VM, so instruction counts come from valgrind/cachegrind. Landau-Vishkin uses AVX-512 where available (scalar fallback otherwise). Not implemented: the longest-common-extension structure for worst-case O(nk) Landau-Vishkin, multithreading, FM-index search, reverse-complement mapping.

## Reproduce

```
make test      # build (g++ -O3 -march=native -std=c++17) and run every exactness test against Sellers; a few minutes
make bench     # the complete pipeline from scratch: build, data, tests, calibration, all benchmarks, plots; about 2 hours on one Xeon core
make quick     # smoke run of the pipeline, a few minutes, tiny samples
scripts/map_reads GENOME.fa READS.fa -k 8 --format sam    # map reads with CIGAR output (both files FASTA)
```

Expected end of `make test` (each script prints one summary line; all must say 0 failures): `300 cases (1398 edlib pattern checks), 0 failures`; `1440 comparisons on real data, 0 failures`; `320 comparisons ... 0 failures` (block-pruned Myers); `600 comparisons ... 0 failures` (Landau-Vishkin); `1120 comparisons ... 0 failures` (generalised filters, 792 reads with a true match); `484805 (split, edit script) pairs checked ... 0 failures` (pigeonhole bound, exhaustive); `3200 random reads in 320 runs ... 0 failures` (alignment output).

Datasets (provenance, URLs, lengths and SHA-256 in `data/SOURCES.txt`): *E. coli* 536 (main benchmark text), *B. fragilis* (blind evaluation only: no constant was fitted on it), SARS-CoV-2, phage lambda, and human and orangutan mitochondrial genomes (real homologs, about 15% apart). Reads are simulated from these texts with `scripts/simulate.py`.

Hardware note: all reported numbers come from one Intel Xeon (2.10 GHz, virtual machine) with AVX-512, and the binary is built with `-march=native`. The Landau-Vishkin (`lv`) speed-ups depend on AVX-512 (a machine without it uses the slower scalar path); Myers, the filters and `filter_auto` do not. The claims in the report are **ratios** (for example speedup over a Myers scan, from interleaved per-read timing), not absolute times, which drift by tens of percent on a shared VM; a different machine will give different absolute times and somewhat different Landau-Vishkin ratios.

Folder layout: `src/` (`approx.cpp`, the single C++ source of all algorithms), `scripts/` (data fetch, simulation, benchmarks, plots, `run_all.sh`), `tests/` (exactness and cross-check tests), `results/` (all tables, CSV and plots produced by the pipeline), `report/` (`REPORT.md`), `data/` (genomes, provenance, temporary files; large files are not tracked by git).
