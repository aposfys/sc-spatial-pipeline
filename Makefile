.PHONY: install data grid analysis test clean clean-data all

PYTHON ?= python3
DATASET ?= visium_hne

all: analysis

## Install the package, the scverse stack and dev tooling. CI installs only ".[dev]",
## which is enough for the tests.
install:
	$(PYTHON) -m pip install -e ".[dev,sc,spatial]"

## Fetch one public spatial dataset, cached
data:
	$(PYTHON) -m scspatial.cli --data-dir data fetch --dataset $(DATASET)

## Run every configuration in the grid, storing labels per configuration
grid: data
	$(PYTHON) -m scspatial.cli --data-dir data grid --dataset $(DATASET)

## Stability table, findings.json and runs.json, all committed
analysis: grid
	$(PYTHON) -m scspatial.cli sensitivity

test:
	$(PYTHON) -m pytest -q

## Removes the grid intermediate only. The committed results are kept.
clean:
	rm -f results/grid.pkl
	find . -name __pycache__ -type d -exec rm -rf {} +

clean-data: clean
	rm -f data/*.h5ad data/*.zarr
