.PHONY: install test lint demo schemas safety clean

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

schemas:
	python scripts/export_schemas.py

safety:
	python scripts/check_repository_safety.py

clean:
	python -c "from pathlib import Path; import shutil; [shutil.rmtree(p, ignore_errors=True) for p in map(Path, ['results', 'build', 'dist', '.pytest_cache', '.mypy_cache', '.ruff_cache'])]"
