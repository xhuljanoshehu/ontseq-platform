"""Structural rules of the package, checked on every run instead of remembered.

Each rule below records a boundary that is true today and that a past defect crossed. The
stage graph used to be re-wired per process by command-line registration, which made the
Desktop service ignore a configured coverage policy and let a stale copy-number report
enter a result. These tests keep that shape from returning quietly: they parse the source,
so they need no tool, no data and no network.

Adding an exception is allowed, but only by editing an allowlist here, in review, with the
reason next to it.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

import pytest

from ontseq_platform.pipeline import stages

PACKAGE_ROOT = Path(__file__).resolve().parents[1] / "src" / "ontseq_platform"
PACKAGE = "ontseq_platform"

#: Modules that compose a run from configuration: they may import the orchestrator.
COMPOSITION_ROOTS = frozenset(
    {
        "ontseq_platform.profile_analysis",
        "ontseq_platform.runtime_cli",
        "ontseq_platform.service.app",
        "ontseq_platform.service.befund",
        "ontseq_platform.system_smoke",
        "ontseq_platform.watchfolder",
    }
)

#: Declarative or storage modules that document themselves as dependency-free.
DEPENDENCY_FREE = frozenset(
    {"ontseq_platform.pipeline.stages", "ontseq_platform.pipeline.envelope"}
)

#: The contract module may import only small, contract-level helpers.
MODELS_ALLOWED_IMPORTS = frozenset({"ontseq_platform.iscn_syntax"})

#: Entry points that may import the command-line modules.
CLI_IMPORTERS = frozenset(
    {
        "ontseq_platform.__main__",
        "ontseq_platform.cnv.__main__",
        "ontseq_platform.entrypoint",
        "ontseq_platform.gatk_adapter.__main__",
    }
)
CLI_MODULES = frozenset(
    {"ontseq_platform.runtime_cli", "ontseq_platform.cli", "ontseq_platform.entrypoint"}
)


@dataclass(frozen=True)
class _Import:
    source: str
    target: str
    #: Imported when the module loads, as opposed to inside a function or TYPE_CHECKING.
    at_import_time: bool
    line: int


def _module_name(path: Path) -> str:
    parts = list(path.relative_to(PACKAGE_ROOT.parent).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _is_type_checking(node: ast.If) -> bool:
    test = node.test
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING"
    )


def _imports_of(path: Path) -> list[_Import]:
    module = _module_name(path)
    package = module if path.name == "__init__.py" else module.rsplit(".", 1)[0]
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found: list[_Import] = []

    def resolve(node: ast.ImportFrom) -> str | None:
        if node.level == 0:
            return node.module
        base = package.split(".")
        base = base[: len(base) - (node.level - 1)]
        return ".".join([*base, node.module] if node.module else base)

    def visit(nodes: list[ast.stmt], *, import_time: bool) -> None:
        for node in nodes:
            if isinstance(node, ast.ImportFrom):
                target = resolve(node)
                if target is None:
                    continue
                names = [target] if node.module else [f"{target}.{a.name}" for a in node.names]
                for name in names:
                    found.append(_Import(module, name, import_time, node.lineno))
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    found.append(_Import(module, alias.name, import_time, node.lineno))
            elif isinstance(node, ast.If):
                guarded = _is_type_checking(node)
                visit(node.body, import_time=import_time and not guarded)
                visit(node.orelse, import_time=import_time)
            elif isinstance(node, ast.Try):
                for block in (node.body, node.orelse, node.finalbody):
                    visit(block, import_time=import_time)
                for handler in node.handlers:
                    visit(handler.body, import_time=import_time)
            elif isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
                visit(
                    node.body,
                    import_time=False if not isinstance(node, ast.ClassDef) else import_time,
                )
            else:
                for child in ast.iter_child_nodes(node):
                    if isinstance(child, ast.stmt):
                        visit([child], import_time=import_time)

    visit(tree.body, import_time=True)
    return found


def _package_modules() -> dict[str, Path]:
    return {
        _module_name(path): path
        for path in sorted(PACKAGE_ROOT.rglob("*.py"))
        if "__pycache__" not in path.parts
    }


MODULES = _package_modules()
IMPORTS = [item for path in MODULES.values() for item in _imports_of(path)]


def _internal(target: str) -> bool:
    return target == PACKAGE or target.startswith(f"{PACKAGE}.")


def _is(target: str, module: str) -> bool:
    return target == module or target.startswith(f"{module}.")


def test_the_parser_sees_the_whole_package() -> None:
    assert "ontseq_platform.pipeline.runner" in MODULES
    assert "ontseq_platform.cnv.lane" in MODULES
    assert any(item.target == "ontseq_platform.pipeline.context" for item in IMPORTS)


def test_only_composition_roots_import_the_orchestrator_when_they_load() -> None:
    """Lanes depend on the stage contract (``pipeline.context``), not on the runner.

    A lane importing the runner at import time is how an orchestrator ends up importing
    itself through every lane it composes.
    """
    offenders = sorted(
        f"{item.source}:{item.line}"
        for item in IMPORTS
        if item.at_import_time
        and _is(item.target, "ontseq_platform.pipeline.runner")
        and item.source not in COMPOSITION_ROOTS
        and item.source != "ontseq_platform.pipeline.runner"
    )
    assert offenders == []


@pytest.mark.parametrize("module", sorted(DEPENDENCY_FREE))
def test_declarative_pipeline_modules_stay_dependency_free(module: str) -> None:
    offenders = sorted(
        item.target
        for item in IMPORTS
        if item.source == module and (_internal(item.target) or _is(item.target, "pydantic"))
    )
    assert offenders == []


def test_the_contract_models_do_not_reach_into_adapters() -> None:
    offenders = sorted(
        item.target
        for item in IMPORTS
        if item.source == "ontseq_platform.models"
        and _internal(item.target)
        and not any(_is(item.target, allowed) for allowed in MODELS_ALLOWED_IMPORTS)
    )
    assert offenders == []


def test_renderers_do_not_depend_on_execution_or_the_service() -> None:
    """Reports render validated contracts; they never drive or query a run."""
    renderers = {
        name
        for name in MODULES
        if name.rsplit(".", 1)[-1].startswith(("report", "workbook")) and name.count(".") == 1
    }
    assert {"ontseq_platform.report", "ontseq_platform.workbook"} <= renderers
    forbidden = ("ontseq_platform.pipeline", "ontseq_platform.service", *CLI_MODULES)
    offenders = sorted(
        f"{item.source} -> {item.target}"
        for item in IMPORTS
        if item.source in renderers and any(_is(item.target, bad) for bad in forbidden)
    )
    assert offenders == []


def test_only_the_runtime_cli_starts_the_service() -> None:
    offenders = sorted(
        f"{item.source} -> {item.target}"
        for item in IMPORTS
        if _is(item.target, "ontseq_platform.service")
        and not item.source.startswith("ontseq_platform.service")
        and item.source != "ontseq_platform.runtime_cli"
    )
    assert offenders == []


def test_only_entry_points_import_command_line_modules() -> None:
    offenders = sorted(
        f"{item.source} -> {item.target}"
        for item in IMPORTS
        if any(item.target == cli for cli in CLI_MODULES)
        and item.source not in CLI_IMPORTERS
        and item.source not in CLI_MODULES
    )
    assert offenders == []


def test_production_code_never_imports_tests() -> None:
    offenders = sorted(
        f"{item.source} -> {item.target}"
        for item in IMPORTS
        if item.target == "tests" or item.target.startswith(("tests.", "test_"))
    )
    assert offenders == []


class _GraphMutation(ast.NodeVisitor):
    """Finds writes to the stage graph: item assignment, deletion or mutating calls."""

    NAMES = frozenset({"IMPLEMENTATIONS", "SPEC_BY_STAGE"})
    MUTATORS = frozenset({"update", "pop", "popitem", "clear", "setdefault", "__setitem__"})

    def __init__(self) -> None:
        self.found: list[int] = []

    def _names_graph(self, node: ast.AST) -> bool:
        if isinstance(node, ast.Name):
            return node.id in self.NAMES
        if isinstance(node, ast.Attribute):
            return node.attr in self.NAMES
        if isinstance(node, ast.Call):
            # ``cast(MutableMapping[...], SPEC_BY_STAGE)`` is how the old registration wrote.
            return any(self._names_graph(argument) for argument in node.args)
        return False

    def _check_target(self, target: ast.AST) -> None:
        if isinstance(target, ast.Subscript) and self._names_graph(target.value):
            self.found.append(target.lineno)

    def visit_Assign(self, node: ast.Assign) -> None:
        for target in node.targets:
            self._check_target(target)
        self.generic_visit(node)

    def visit_AugAssign(self, node: ast.AugAssign) -> None:
        self._check_target(node.target)
        self.generic_visit(node)

    def visit_Delete(self, node: ast.Delete) -> None:
        for target in node.targets:
            self._check_target(target)
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        function = node.func
        if (
            isinstance(function, ast.Attribute)
            and function.attr in self.MUTATORS
            and self._names_graph(function.value)
        ):
            self.found.append(node.lineno)
        self.generic_visit(node)


def test_no_module_rewires_the_stage_graph() -> None:
    """A run-specific difference is configuration, never a process-wide replacement."""
    offenders: list[str] = []
    for name, path in MODULES.items():
        visitor = _GraphMutation()
        visitor.visit(ast.parse(path.read_text(encoding="utf-8")))
        offenders.extend(f"{name}:{line}" for line in visitor.found)
    assert offenders == []


def test_the_stage_specifications_are_read_only_at_run_time() -> None:
    assert isinstance(stages.SPEC_BY_STAGE, MappingProxyType)
    with pytest.raises(TypeError):
        stages.SPEC_BY_STAGE[stages.StageId.CNV] = stages.SPEC_BY_STAGE[  # type: ignore[index]
            stages.StageId.QC
        ]


def test_no_production_module_keeps_process_global_state_by_rebinding() -> None:
    """``global`` is how per-process lane settings used to be installed.

    Mutable module state makes a long-running service depend on which command ran first.
    """
    offenders = sorted(
        f"{name}:{node.lineno}"
        for name, path in MODULES.items()
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8")))
        if isinstance(node, ast.Global)
    )
    assert offenders == []
