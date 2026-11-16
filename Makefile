.PHONY: bench quick test clean
bench:
	./scripts/run_all.sh
quick:
	./scripts/run_all.sh --quick
test:
	g++ -O3 -march=native -std=c++17 -o approx src/approx.cpp && python3 tests/edlib_semantics.py && python3 tests/crosscheck.py
clean:
	rm -f approx approx_generic results/*.png results/*.csv results/summary.md data/_*.txt
