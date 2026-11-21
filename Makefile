.PHONY: bench quick test clean
bench:
	./scripts/run_all.sh
quick:
	./scripts/run_all.sh --quick
test:
	g++ -O3 -march=native -std=c++17 -o approx src/approx.cpp && python3 tests/edlib_semantics.py && python3 tests/crosscheck.py && python3 tests/crosscheck_real.py && python3 tests/test_pruned.py && python3 tests/test_lv.py && python3 tests/test_filter_s.py && python3 tests/test_mapping.py
clean:
	rm -f approx approx_generic results/*.png results/*.csv results/*.md results/*.json data/_*.txt
