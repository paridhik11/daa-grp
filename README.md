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
| `src/approx.cpp` | Sellers, banded, Myers, Myers with block pruning, Landau-Vishkin (direct extension), equal-split filter, frequency-aware filter (hits, and hits + lambda * windows), BLAST-style heuristic, traceback (`map`) |
| `scripts/` | `fetch_data.py`, `simulate.py`, `calibrate.py`, `run_bench.py`, `plot.py`, `instr_profile.py`, `expected_hits.py`, `real_divergence.py`, `objective_test.py`, `memory.py`, `fig_dp_matrix.py`, `mapping_example.py`, `map_reads`, `run_all.sh` |
| `tests/` | `crosscheck.py`, `crosscheck_real.py`, `test_pruned.py`, `test_mapping.py` (all against Sellers), `edlib_semantics.py` |
| `data/SOURCES.txt` | Dataset provenance (URL, header, length, GC, SHA-256) |

Datasets: E. coli 536 (main), B. fragilis (blind evaluation only), SARS-CoV-2, phage lambda, human and orangutan mtDNA, fetched from GitHub raw files. K-12 MG1655 could not be obtained (NCBI, ENA, UCSC and Ensembl were unreachable), so no dataset is called K-12. Edlib and WFA2 return the best match per read, not all matches within k. Hardware performance counters are unavailable in the sandbox VM, so instruction counts come from valgrind/cachegrind. Not implemented: the longest-common-extension structure for worst-case O(nk) Landau-Vishkin, multithreading, FM-index search, reverse-complement mapping.
