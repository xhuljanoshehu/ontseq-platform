.PHONY: install test lint demo local-smoke schemas safety clean

install:
	python -m pip install -e ".[dev]"

test:
	python -m unittest discover -s tests -v

lint:
	ruff check .
	ruff format --check .
	mypy src

demo:
	python -m ontseq_platform demo --output-dir results/demo

local-smoke:
	PYTHONPATH=src python -m ontseq_platform local-smoke --output-dir results/local-smoke

schemas:
	python scripts/export_schemas.py

safety:
	python scripts/check_repository_safety.py

clean:
	python -c "from pathlib import Path; import shutil; [shutil.rmtree(p, ignore_errors=True) for p in map(Path, ['results', 'build', 'dist', '.pytest_cache', '.mypy_cache', '.ruff_cache'])]"
