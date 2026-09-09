PY ?= python3
.PHONY: all quick test figures clean
all:            ## every number and figure (about 40 min)
	$(PY) run_all.py
quick:          ## everything except the four long searches
	$(PY) run_all.py --quick
test:           ## compare results/*.json with the values quoted in the manuscript
	$(PY) -m pytest -q tests
figures:        ## only the eight figures (requires results of the analyses)
	cd scripts && for f in fig1_invariance fig2_sites fig3_plateau fig4_plane fig5_integer fig6_uniqueness fig7_tests fig8_budget; do $(PY) $$f.py; done
clean:
	rm -rf results/*.json results/*.log results/*.npy results/*.npz figures/*.pdf S1_Data/*.json S1_Data/*.csv esbm/__pycache__ scripts/__pycache__ tests/__pycache__ .pytest_cache
