"""Curated research capabilities, NOT an executable provider registry.

Existing runtime adapters are left unchanged. A catalog entry does not establish
analytical validity, installation, or compatibility with a particular specimen.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

CATALOG_VERSION = "0.1.0"


@dataclass(frozen=True)
class CallerDefinition:
    caller_id: str
    label: str
    analyses: frozenset[str]
    signal_families: frozenset[str]
    method_family: str
    existing_adapter: bool
    source_repository: str


_CALLERS = (
    CallerDefinition(
        "qdnaseq_ace",
        "QDNAseq + ACE",
        frozenset({"cnv", "purity_ploidy"}),
        frozenset({"depth"}),
        "qdnaseq_ace",
        True,
        "https://github.com/tgac-vumc/ACE",
    ),
    CallerDefinition(
        "ont_spectre",
        "ONT-Spectre",
        frozenset({"cnv"}),
        frozenset({"depth", "snv"}),
        "spectre",
        False,
        "https://github.com/nanoporetech/ont-spectre",
    ),
    CallerDefinition(
        "ichorcna",
        "ichorCNA",
        frozenset({"cnv", "purity_ploidy"}),
        frozenset({"depth"}),
        "ichorcna",
        False,
        "https://github.com/broadinstitute/ichorCNA",
    ),
    CallerDefinition(
        "sniffles2",
        "Sniffles2",
        frozenset({"sv"}),
        frozenset({"breakpoint"}),
        "sniffles",
        True,
        "https://github.com/fritzsedlazeck/Sniffles",
    ),
    CallerDefinition(
        "cutesv",
        "cuteSV",
        frozenset({"sv"}),
        frozenset({"breakpoint"}),
        "cutesv",
        True,
        "https://github.com/tjiangHIT/cuteSV",
    ),
    CallerDefinition(
        "savana",
        "SAVANA",
        frozenset({"sv", "cnv", "allelic_cn", "purity_ploidy"}),
        frozenset({"breakpoint", "depth", "snv"}),
        "savana",
        False,
        "https://github.com/cortes-ciriano-lab/savana",
    ),
    CallerDefinition(
        "severus",
        "Severus",
        frozenset({"sv", "complex_sv"}),
        frozenset({"breakpoint", "phasing"}),
        "severus",
        False,
        "https://github.com/KolmogorovLab/Severus",
    ),
    CallerDefinition(
        "wakhan",
        "Wakhan",
        frozenset({"cnv", "allelic_cn", "purity_ploidy"}),
        frozenset({"depth", "phasing", "breakpoint"}),
        "wakhan",
        False,
        "https://github.com/KolmogorovLab/Wakhan",
    ),
)

CALLER_CATALOG: Mapping[str, CallerDefinition] = MappingProxyType(
    {entry.caller_id: entry for entry in _CALLERS}
)
