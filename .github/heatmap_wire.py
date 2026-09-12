"""One-shot call-site patches for the methylation heatmap. Deleted after success."""

from pathlib import Path


def patch(path: str, old: str, new: str) -> None:
    target = Path(path)
    text = target.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected 1 match, found {count}")
    target.write_text(text.replace(old, new, 1), encoding="utf-8")
    print(f"patched {path}")


patch(
    "src/ontseq_platform/pipeline/runner.py",
    """    render_html(
        result,
        ctx.envelope.path(ctx.path(REPORT_HTML)),
        target_coverage=target_coverage,
        selection_coverage=selection_coverage,
        qc_histogram=qc_histogram,
    )
""",
    """    render_html(
        result,
        ctx.envelope.path(ctx.path(REPORT_HTML)),
        target_coverage=target_coverage,
        selection_coverage=selection_coverage,
        qc_histogram=qc_histogram,
        methylation_report=load_methylation_report(ctx),
    )
""",
)

patch(
    "src/ontseq_platform/cnv/extension.py",
    """from ..pipeline.stages import SPEC_BY_STAGE, StageId, StageSpec, VerificationStatus
from ..report import render_html
from ..report_plots import CnvChromosomeBar, cnv_genome_svg
""",
    """from ..pipeline.stages import SPEC_BY_STAGE, StageId, StageSpec, VerificationStatus
from ..qc import read_length_histogram_from_tsv
from ..report import render_html
from ..report_plots import CnvChromosomeBar, ReadLengthBin, cnv_genome_svg
""",
)

patch(
    "src/ontseq_platform/cnv/extension.py",
    """    render_html(
        result,
        html_path,
        target_coverage=target_coverage,
        selection_coverage=selection_coverage,
    )
""",
    """    histogram_path = ctx.envelope.path(pipeline_runner.QC_READ_LENGTH_HISTOGRAM)
    qc_histogram = (
        [
            ReadLengthBin(start=start, end=end, count=count, bases=bases)
            for start, end, count, bases in read_length_histogram_from_tsv(
                histogram_path.read_text(encoding="utf-8")
            )
        ]
        if histogram_path.is_file()
        else None
    )
    render_html(
        result,
        html_path,
        target_coverage=target_coverage,
        selection_coverage=selection_coverage,
        qc_histogram=qc_histogram,
        methylation_report=pipeline_runner.load_methylation_report(ctx),
    )
""",
)

patch(
    "CHANGELOG.md",
    """  Pure presentation of the normalized QDNAseq/ACE report; no fit, threshold, contract
  or reportability change.

### Fixed
""",
    """  Pure presentation of the normalized QDNAseq/ACE report; no fit, threshold, contract
  or reportability change.
- Render the normalized methylation report as a deterministic inline-SVG heatmap:
  regions in genome order, one row per modification code, measured fractions on a
  0-1 blue scale, and coverage-floor misses as hatched not-measurable cells rather
  than zeros. Failed and no-call totals from summary metrics stay visible beside
  the figure. Presentation only: no pileup, threshold, contract or reportability
  change.

### Fixed
""",
)
