.PHONY: install data grid analysis test clean clean-data all

PYTHON ?= python3
DATASET ?= visium_hne

all: analysis

## Install the package plus dev tooling. The scverse stack is an optional extra:
##   pip install -e ".[sc,spatial]"
install:
	$(PYTHON) -m pip install -e ".[dev]"

## Fetch one public spatial dataset, cached
data:
	$(PYTHON) -m scspatial.cli fetch --dataset $(DATASET)

## Run every configuration in the grid, storing labels per configuration
grid: data
	$(PYTHON) -m scspatial.cli grid --dataset $(DATASET)

## Stability across configurations, including rare-population retention
analysis: grid
	$(PYTHON) -m scspatial.cli sensitivity

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf results/*
	find . -name __pycache__ -type d -exec rm -rf {} +

clean-data: clean
	rm -f data/*.h5ad data/*.zarr
