# Approximate DNA matching: comparative study (DAA project, Problem 8)

One command (build, fetch real genomes, correctness tests, NNLS cost-model calibration, instruction profile, benchmarks incl. a 100-read headline experiment with median and 10th-90th percentile, stress tests, plots, summary; about 35 minutes):

```
make bench
```

Smoke run (about 2 minutes): `make quick`. Correctness tests only: `make test`.

| Path | Contents |
|---|---|
| `report/REPORT.md` | Full report: problem, algorithms, proofs, complexity, conclusion, gaps, references |
| `results/summary.md` | All tables (break-even, cost model, cells per second, stress, datasets) |
| `results/*.png` | Plots; `results/bench.csv` raw data |
| `src/approx.cpp` | Sellers, banded, Myers, equal-split filter, frequency-aware filter, BLAST-style heuristic |
| `scripts/` | `fetch_data.py`, `simulate.py`, `calibrate.py`, `instr_profile.py`, `run_bench.py`, `plot.py`, `run_all.sh` |
| `tests/` | `crosscheck.py` (vs Sellers and Edlib), `edlib_semantics.py` |
| `data/SOURCES.txt` | Dataset provenance (URL, header, length, GC, SHA-256) |

Datasets: E. coli 536 (main), SARS-CoV-2, phage lambda, human mtDNA, fetched from GitHub raw files. K-12 MG1655 could not be obtained (NCBI, ENA, UCSC and Ensembl were unreachable), so no dataset is called K-12. Edlib and WFA2 return the best match per read, not all matches within k. Hardware performance counters are unavailable in the sandbox VM, so instruction counts come from valgrind/cachegrind.
