# Convenience targets. Every one of these is a single plain command -- see the
# README if you would rather run them directly.

.PHONY: help install demo demo-synthetic figure test all clean

help:
	@echo "make install  - install dependencies"
	@echo "make demo     - analyse the measured spectra in data/, print the tables"
	@echo "make demo-synthetic - same, on the synthetic pair with known answers"
	@echo "make figure   - render the annotated figure used in the README"
	@echo "make test     - run the ground-truth eval harness"
	@echo "make all      - demo + figure + test"

install:
	pip install -r requirements.txt

demo:
	python3 analyse_raman.py --outdir results

demo-synthetic:
	python3 analyse_raman.py examples/synthetic/*.txt --outdir /tmp/raman-synthetic

figure:
	python3 plot_spectrum.py --out results/spectrum.png

test:
	python3 -m pytest tests/ -q

all: demo figure test

# Note: results/ is committed on purpose -- the README embeds results/spectrum.png
# and links the tables, so a reviewer sees output before running anything.
# `make all` regenerates it; clean only removes caches.
clean:
	rm -rf __pycache__ tests/__pycache__ .pytest_cache
