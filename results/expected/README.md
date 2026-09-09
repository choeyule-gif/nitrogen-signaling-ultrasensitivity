# results/expected

The `results/*.json` files as produced on the machine the manuscript was written from
(CPython 3.14.6, numpy 2.5.2, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0; macOS 26, Apple silicon).
`python run_all.py` regenerates `results/*.json`; `tests/test_reproduce.py` compares the fresh
files with the manuscript's quoted values, not with these copies.  These copies let a reader
check a single number without running anything.
