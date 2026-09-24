"""Presentation of the normalized modified-base report, without adding any interpretation.

Every value shown here already exists in :class:`~ontseq_platform.methylation.MethylationReport`.
Three presentation rules follow from that report's own contract and are not negotiable:

* a region below the coverage floor is shown as *not measurable*, never as ``0``: its
  fraction is ``None`` in the report and stays empty in the workbook;
* failed-threshold, no-call, deletion and diff calls stay visible next to the valid
  denominator, so "not measured" cannot be read as "not modified";
* nothing is rendered for a report that does not belong to the result being presented.
"""

from __future__ import annotations

from collections.abc import Sequence

from .methylation import MODKIT_MODIFICATION_NAMES, MethylationRegionSummary, MethylationReport
from .models import PipelineResult
from .report_plots import chromosome_sort_key

#: Shown wherever a region has no site at the coverage floor.
NOT_MEASURABLE = "not measurable (below the coverage floor)"

_CODE_RANK = {"m": 0, "h": 1, "a": 2}

REGION_HEADERS: tuple[str, ...] = (
    "Region",
    "Chromosome",
    "Start (0-based)",
    "End",
    "Modification",
    "Modification name",
    "Measurement",
    "Sites",
    "Sites at coverage floor",
    "Valid calls",
    "Modified calls",
    "Canonical calls",
    "Other-modification calls",
    "Failed-threshold calls",
    "No-call calls",
    "Deletion calls",
    "Diff calls",
    "Mean modified fraction (call-weighted)",
    "Median site modified fraction",
    "Mean valid coverage",
)


def validate_methylation_identity(result: PipelineResult, report: MethylationReport | None) -> None:
    """Refuse to present a methylation report that belongs to another result.

    The report must describe the same sample and build, be the pileup whose bedMethyl
    checksum the result recorded, and agree with the result's methylation module outcome.
    """
    if report is None:
        return
    if (report.sample_id, report.genome_build) != (
        result.manifest.sample_id,
        result.manifest.assay.genome_build,
    ):
        raise ValueError("methylation report sample/build identity differs from the result")
    if (
        result.provenance.reference_checksums.get("bedmethyl")
        != report.bedmethyl_fingerprint.sha256
    ):
        raise ValueError("methylation report is not the pileup the result recorded")
    if not any(
        item.module.value == "methylation"
        and item.status == report.status
        and report.tool in item.tools
        for item in result.modules
    ):
        raise ValueError("methylation report contradicts the normalized module outcome")


def ordered_regions(report: MethylationReport) -> list[MethylationRegionSummary]:
    """Regions in genome order, then by modification code (5mC, 5hmC, 6mA)."""
    return sorted(
        report.regions,
        key=lambda region: (
            chromosome_sort_key(region.chromosome),
            region.start is None,
            region.start or 0,
            region.end or 0,
            region.region_id,
            _CODE_RANK.get(region.modification_code.value, 50),
        ),
    )


def modification_label(region: MethylationRegionSummary) -> str:
    return MODKIT_MODIFICATION_NAMES[region.modification_code]


def measurement_label(region: MethylationRegionSummary) -> str:
    return "measured" if region.mean_modified_fraction is not None else NOT_MEASURABLE


def region_rows(report: MethylationReport) -> list[list[object]]:
    """One workbook row per region and modification code; ``None`` stays empty."""
    return [
        [
            region.region_id,
            region.chromosome,
            region.start,
            region.end,
            modification_label(region),
            region.modification_name,
            measurement_label(region),
            region.sites_total,
            region.sites_at_minimum_coverage,
            region.valid_call_count,
            region.modified_call_count,
            region.canonical_call_count,
            region.other_mod_call_count,
            region.fail_call_count,
            region.nocall_call_count,
            region.delete_call_count,
            region.diff_call_count,
            region.mean_modified_fraction,
            region.median_site_modified_fraction,
            region.mean_valid_coverage,
        ]
        for region in ordered_regions(report)
    ]


def _tagged_reads(report: MethylationReport) -> str:
    if report.reads_with_modified_base_tags is not None:
        return str(report.reads_with_modified_base_tags)
    if report.policy.verify_modified_base_tags:
        return "unknown: the installed samtools could not evaluate the tag probe"
    return "not probed: disabled by the policy"


def methylation_facts(report: MethylationReport | None) -> list[list[object]]:
    """Field/value rows describing what the pileup measured and under which policy."""
    if report is None:
        return [
            ["Status", "NOT_RECORDED"],
            [
                "Meaning",
                "No current methylation report belongs to this result; see the module "
                "status. This is not an unmethylated result.",
            ],
        ]
    policy = report.policy
    metrics = report.summary_metrics
    measured = sum(region.sites_at_minimum_coverage for region in report.regions)
    return [
        ["Status", report.status.value],
        [
            "Meaning",
            "Descriptive modified-base fractions; no methylation threshold, region set or "
            "classifier here is analytically validated."
            if report.status.value == "COMPLETED"
            else "No site reached the coverage floor: the sample was not measurable at this "
            "depth. This NO_CALL does not say the DNA is unmethylated.",
        ],
        ["Region source", report.region_source.value],
        ["Region rows", len(report.regions)],
        ["Site observations at coverage floor", measured],
        ["Pileup sites (all codes)", metrics.get("site_count")],
        ["Excluded non-canonical contig sites", metrics.get("skipped_non_canonical_site_count")],
        ["Failed-threshold calls (all sites)", metrics.get("fail_call_count")],
        ["No-call calls (all sites)", metrics.get("nocall_call_count")],
        ["Reads with MM tags", _tagged_reads(report)],
        ["Policy", f"{policy.profile_id} ({policy.status})"],
        [
            "Modification codes",
            ", ".join(MODKIT_MODIFICATION_NAMES[c] for c in policy.modification_codes),
        ],
        ["CpG only / strands combined", f"{policy.cpg_only} / {policy.combine_strands}"],
        ["Call confidence threshold (pinned)", policy.filter_threshold],
        ["Minimum valid coverage per site", policy.minimum_valid_coverage],
        ["Tool", f"{report.tool.name} {report.tool.version}"],
        ["bedMethyl SHA-256", report.bedmethyl_fingerprint.sha256],
        [
            "Target BED SHA-256",
            report.target_bed_fingerprint.sha256 if report.target_bed_fingerprint else None,
        ],
        ["Genome build", report.genome_build.value],
        ["Warnings", "\n".join(report.warnings)],
        ["Limitations", "\n".join(report.limitations)],
    ]


def fraction_text(value: float | None) -> str:
    return NOT_MEASURABLE if value is None else f"{value:.1%}"


def coverage_text(value: float | None) -> str:
    return "not available" if value is None else f"{value:.1f}×"


def locus_text(region: MethylationRegionSummary) -> str:
    if region.start is None or region.end is None:
        return f"{region.chromosome} (whole chromosome, observed sites only)"
    return f"{region.chromosome}:{region.start}-{region.end}"


def html_rows(regions: Sequence[MethylationRegionSummary]) -> list[tuple[str, ...]]:
    """Plain-text cells for the HTML table; the renderer escapes every one of them."""
    return [
        (
            region.region_id,
            locus_text(region),
            modification_label(region),
            f"{region.sites_at_minimum_coverage} / {region.sites_total}",
            str(region.valid_call_count),
            str(region.modified_call_count),
            fraction_text(region.mean_modified_fraction),
            fraction_text(region.median_site_modified_fraction),
            coverage_text(region.mean_valid_coverage),
            f"{region.fail_call_count} / {region.nocall_call_count}",
        )
        for region in regions
    ]
