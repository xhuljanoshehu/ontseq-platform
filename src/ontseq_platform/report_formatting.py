"""Presentation-only escaping and conservative, allowlisted tool metadata export."""

from __future__ import annotations

import html
import json
import math
import re
from collections.abc import Mapping

REDACTED = "[redacted path-like token]"
# Do not match a closing HTML tag or biological slash terms (CNV/SV, tumor/normal).
_PATH = re.compile(
    r"(?<![\w<])(?:file://|[A-Za-z]:[\\/]|\\\\|/(?![/>]))[^\s<>\"']+"
    r"|\b(?:[\w.-]+/)+[\w.-]+\.(?:bam|bai|csi|sam|bed|vcf|fa|fasta|fai|"
    r"yaml|yml|json|tsv|txt|pod5|fast5)(?:\.gz)?\b",
    re.IGNORECASE,
)
_PATH_FIELD = re.compile(
    r"[\"']?(?:input_path|output_path|source_path|bam_path|reference_path|model_path)"
    r"[\"']?\s*[:=]\s*(?:\"[^\"]*\"|'[^']*'|[^\s,;}]+)",
    re.IGNORECASE,
)
_QUOTED_PATH = re.compile(r"""(["'])(?:file://|[A-Za-z]:[\\/]|\\\\|/)[^\n]*?\1""")


def redact_paths(value: str) -> str:
    value = _PATH_FIELD.sub(REDACTED, value)
    value = _QUOTED_PATH.sub(REDACTED, value)
    return _PATH.sub(REDACTED, value)


def cell(value: object) -> str:
    return html.escape(redact_paths("" if value is None else str(value)))


def optional_number(value: float | int | None, *, percent: bool = False, places: int = 1) -> str:
    if value is None or not math.isfinite(value):
        return "not available"
    return f"{value * (100 if percent else 1):.{places}f}" + ("%" if percent else "×")


_COMMON = {"threads", "expected_version", "adapter_verification"}
_CNV = {
    "profile",
    "bin_sizes_kbp",
    "bins_kb",
    "primary_bin_size_kbp",
    "ace_penalty",
    "penalty",
    "ploidy_min",
    "ploidy_max",
    "ploidy_step",
}
_TOOL_PARAMETERS = {
    "sniffles2": _COMMON
    | {
        "minsupport",
        "minsvlen",
        "mapq",
        "caller_pass_only",
        "normalizer_pass_only",
        "symbolic",
        "mosaic",
        "output_read_names",
    },
    "cutesv": _COMMON
    | {
        "min_support",
        "min_size",
        "max_cluster_bias_INS",
        "diff_ratio_merging_INS",
        "max_cluster_bias_DEL",
        "diff_ratio_merging_DEL",
        "normalizer_pass_only",
    },
    "mosdepth": _COMMON | {"no_per_base", "thresholds", "mapq", "exclude_flags", "target_bed_role"},
    "qdnaseq": _COMMON | _CNV,
    "ace": _COMMON | _CNV,
    "qdnaseq+ace": _COMMON | _CNV,
    "cramino": _COMMON | {"format", "read_length_histogram_requested"},
    "minimap2": _COMMON | {"preset", "md_tag", "soft_clip_supplementary", "carry_fastq_tags"},
    "samtools": _COMMON | {"checks"},
    "dorado": _COMMON
    | {"model", "model_sha256", "device", "modified_bases", "minimum_qscore", "emit_moves"},
    "modkit": _COMMON
    | {
        "filter_threshold",
        "cpg_only",
        "combine_strands",
        "ignored_codes",
        "ontseq_region_assignment",
        "modification_codes",
    },
}


def _safe_parameter(value: object) -> bool:
    if value is None or isinstance(value, bool | int):
        return True
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, str):
        return redact_paths(value) == value and not any(c in value for c in ("\n", "\r", "\x00"))
    if isinstance(value, list | tuple):
        return all(_safe_parameter(item) for item in value)
    return False


def tool_parameters(name: str, parameters: Mapping[str, object]) -> str:
    allowed = _TOOL_PARAMETERS.get(name.lower(), set())
    exported = {
        key: value for key, value in parameters.items() if key in allowed and _safe_parameter(value)
    }
    text = json.dumps(exported, sort_keys=True, allow_nan=False)
    if len(exported) != len(parameters):
        text += " · unsafe or path-like parameters redacted"
    return text
