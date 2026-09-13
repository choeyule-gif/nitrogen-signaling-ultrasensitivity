PY ?= python3
.PHONY: all quick test figdata s1data clean
all:            ## every number, the figure data and the S1 Data package (about 45 min)
	$(PY) run_all.py
quick:          ## the same without the four long searches
	$(PY) run_all.py --quick
test:           ## compare results/*.json with the values quoted in the manuscript
	$(PY) tests/test_reproduce.py
figdata:        ## only the CSV files the MATLAB figures read (needs results/ from a previous run)
	$(PY) scripts/export_figure_data.py && $(PY) scripts/export_extra_figure_data.py
s1data:         ## assemble S1_Data/ for submission
	$(PY) scripts/make_s1_data.py
clean:
	rm -rf results/*.json results/*.log results/*.npy results/*.npz results/panels figures S1_Data \
	       esbm/__pycache__ scripts/__pycache__ tests/__pycache__ .pytest_cache

.PHONY: matched-models
matched-models:
	$(PY) extras/matched_models/run.py

.PHONY: finite-mechanisms
finite-mechanisms:
	$(PY) extras/finite_mechanism/run.py

.PHONY: regulatory-states
regulatory-states:
	$(PY) extras/regulatory_states/code/analyze.py

.PHONY: robust-design
robust-design:
	$(PY) extras/robust_design/code/analyze.py

.PHONY: closed-cascade
closed-cascade:
	$(PY) extras/closed_cascade/code/analyze.py
	$(PY) extras/closed_cascade/code/topology.py
	$(PY) extras/closed_cascade/code/total_input.py
	$(PY) extras/closed_cascade/code/verify_algebra.py
