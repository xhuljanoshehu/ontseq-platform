"""Which guideline criteria the panel design can reach at all, before any caller is wired in.

A criteria table says what a guideline asks. A panel BED says what the sequencer was told to
enrich. Nothing in this repository compared the two, so a gap could sit between them
indefinitely: a criterion can be perfectly implemented, run without error, and still be
answering a question about a gene the run never sequenced deeply.

That gap is real in the shipped design. ``aml_fusion_adaptive_sampling.grch38.buffered.bed``
is a fusion panel. It targets NPM1, FLT3 and RUNX1, and it does not target CEBPA, TP53, or
eight of the nine myelodysplasia-related genes. On the laboratory's measured run the panel
carried roughly 80x and everything outside it roughly 8.8x, so "not in the panel" is not a
detail of bookkeeping -- it is the difference between a depth that supports somatic
small-variant calling and one that does not.

Scope, deliberately narrow. This module answers one question: *for criteria that depend on
small-variant calling, is every gene the criterion names inside the panel?* It does not
second-guess criteria detectable by copy number or by breakpoint. Copy number is read from
the off-target background by design, and the panel is where breakpoints are sought, so panel
membership is not the gating question for those and this module says so rather than
inventing a verdict.

The verdict that matters most is :data:`PARTIAL`. A criterion naming nine genes with one of
them targeted is the dangerous case, not a partial success: evaluated against the genes that
happen to be present it can report "no myelodysplasia-related mutation" having never looked
at eight of them. That is the same bias :func:`ontseq_platform.guideline_criteria.
risk_group_determinable` refuses elsewhere -- absent evidence read as evidence of absence,
skewed towards favourable -- so :data:`PARTIAL` is reported as unreachable, loudly, and never
folded in with the reachable ones.

Nothing here validates a panel, a threshold or a gene identity, and no coordinate is checked
against a gene model. Research use only.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ontseq_platform.guideline_criteria import Criterion

SMALL_VARIANT = "small_variant"

REACHABLE = "all_named_genes_in_panel"
PARTIAL = "some_named_genes_missing_from_panel"
UNREACHABLE = "no_named_gene_in_panel"
NO_GENES_NAMED = "criterion_names_no_genes"
NOT_PANEL_GATED = "not_gated_by_panel_membership"


class PanelReachabilityError(ValueError):
    """Raised when a panel BED cannot be read as a set of named targets."""


@dataclass(frozen=True)
class PanelTarget:
    """One BED interval. Half-open, as BED defines it."""

    chrom: str
    start: int
    end: int
    label: str

    @property
    def span(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Panel:
    path: Path
    targets: tuple[PanelTarget, ...]

    @property
    def labels(self) -> tuple[str, ...]:
        """Every target label in file order, duplicates kept.

        Deliberately parallel to :func:`ontseq_platform.pipeline.panel_lock.target_labels`.
        That module is the authority on panel labels, but it imports yaml to read the lock
        beside the BED, and this module is kept dependency-free so it can run wherever the
        criteria table can. A test pins the two parsers to the same answer.
        """
        return tuple(target.label for target in self.targets)

    @property
    def gene_names(self) -> frozenset[str]:
        return frozenset(self.labels)

    @property
    def total_span(self) -> int:
        return sum(target.span for target in self.targets)

    def targets_for(self, gene: str) -> tuple[PanelTarget, ...]:
        return tuple(target for target in self.targets if target.label == gene)


def load_panel(path: Path) -> Panel:
    """Read a 4-column BED. Comments and blank lines are skipped; a bad row is an error.

    An unnamed or non-numeric interval is refused rather than dropped. A silently shortened
    panel would understate exactly what this module exists to measure.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        raise PanelReachabilityError(f"panel BED is unreadable: {path}") from error

    targets: list[PanelTarget] = []
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 4 or not fields[3].strip():
            raise PanelReachabilityError(f"{path.name}: line {number} has no target label")
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError as error:
            raise PanelReachabilityError(
                f"{path.name}: line {number} has non-numeric coordinates"
            ) from error
        if end <= start:
            raise PanelReachabilityError(
                f"{path.name}: line {number} is empty or reversed ({start}, {end})"
            )
        targets.append(PanelTarget(chrom=fields[0], start=start, end=end, label=fields[3].strip()))

    if not targets:
        raise PanelReachabilityError(f"panel BED holds no intervals: {path}")
    return Panel(path=path, targets=tuple(targets))


@dataclass(frozen=True)
class CriterionReach:
    """Whether the panel puts every gene a criterion names under enrichment."""

    record_id: str
    display_name: str
    status: str
    needs_small_variants: bool
    genes_in_panel: tuple[str, ...]
    genes_missing: tuple[str, ...]

    @property
    def reachable(self) -> bool:
        """True only when the panel can support the criterion as written.

        :data:`PARTIAL` is false. A criterion answered from a subset of the genes it names is
        not a partial answer, it is a wrong one biased towards favourable.
        """
        return self.status in (REACHABLE, NOT_PANEL_GATED)

    def reason(self) -> str:
        if self.status == NOT_PANEL_GATED:
            return (
                "Detectable by copy number or breakpoint, which panel membership does not "
                "gate. This module states no verdict."
            )
        if self.status == REACHABLE:
            return f"Every named gene is a panel target: {', '.join(self.genes_in_panel)}."
        if self.status == NO_GENES_NAMED:
            return "Needs small variants but names no gene, so reachability cannot be decided."
        if self.status == UNREACHABLE:
            return (
                f"No named gene is a panel target ({', '.join(self.genes_missing)}); at "
                "off-target depth this criterion cannot be evaluated."
            )
        return (
            f"{len(self.genes_missing)} of "
            f"{len(self.genes_in_panel) + len(self.genes_missing)} named genes are not panel "
            f"targets ({', '.join(self.genes_missing)}). Evaluating the criterion from the "
            f"remaining {', '.join(self.genes_in_panel)} would report absence for genes that "
            "were never sequenced deeply enough to see."
        )


def assess_criterion(criterion: Criterion, panel: Panel) -> CriterionReach:
    """Decide one criterion against one panel, without consulting its ``assay_status``.

    The criteria table's own ``assay_status`` records that no small-variant caller is wired
    in. That is a statement about the software. This is a statement about the design, and the
    two are worth keeping apart: wiring a caller in changes the first and not the second.
    """
    needs = SMALL_VARIANT in criterion.detectable_by
    present = tuple(gene for gene in criterion.genes if gene in panel.gene_names)
    missing = tuple(gene for gene in criterion.genes if gene not in panel.gene_names)

    if not needs:
        status = NOT_PANEL_GATED
    elif not criterion.genes:
        status = NO_GENES_NAMED
    elif not missing:
        status = REACHABLE
    elif not present:
        status = UNREACHABLE
    else:
        status = PARTIAL

    return CriterionReach(
        record_id=criterion.record_id,
        display_name=criterion.display_name,
        status=status,
        needs_small_variants=needs,
        genes_in_panel=present,
        genes_missing=missing,
    )


@dataclass(frozen=True)
class ReachabilityReport:
    panel_path: Path
    bundle_id: str
    reaches: tuple[CriterionReach, ...]

    def by_status(self, status: str) -> tuple[CriterionReach, ...]:
        return tuple(item for item in self.reaches if item.status == status)

    @property
    def small_variant_criteria(self) -> tuple[CriterionReach, ...]:
        return tuple(item for item in self.reaches if item.needs_small_variants)

    @property
    def blocked(self) -> tuple[CriterionReach, ...]:
        """Small-variant criteria the panel design cannot support, worst case first."""
        order = {PARTIAL: 0, UNREACHABLE: 1, NO_GENES_NAMED: 2}
        blocked = [item for item in self.small_variant_criteria if not item.reachable]
        return tuple(sorted(blocked, key=lambda item: (order[item.status], item.record_id)))

    @property
    def genes_to_add(self) -> tuple[str, ...]:
        """Every gene a small-variant criterion names that the panel does not target.

        This is the list a laboratory needs in order to decide whether to extend the design.
        Adaptive sampling costs no reagent per added region, so the decision is about depth
        dilution across a larger territory, not about consumables -- but it is an assay design
        decision and nothing here makes it.
        """
        missing: set[str] = set()
        for item in self.small_variant_criteria:
            missing.update(item.genes_missing)
        return tuple(sorted(missing))


def assess_panel(
    criteria: Sequence[Criterion], panel: Panel, *, bundle_id: str
) -> ReachabilityReport:
    return ReachabilityReport(
        panel_path=panel.path,
        bundle_id=bundle_id,
        reaches=tuple(assess_criterion(criterion, panel) for criterion in criteria),
    )


def _kb(span: int) -> str:
    return f"{span / 1000:.1f} kb"


def _mb(span: int) -> str:
    return f"{span / 1_000_000:.2f} Mb"


def format_report(report: ReachabilityReport, panel: Panel) -> str:
    """Render the assessment as Markdown, generated so it cannot drift from the inputs."""
    lines = [
        "# Panel reachability of the guideline criteria",
        "",
        "GENERATED FILE - do not edit by hand. Regenerate with:",
        "",
        "    python -m ontseq_platform.panel_reachability",
        "",
        f"Panel: `{report.panel_path.name}` - {len(panel.targets)} targets, "
        f"{_mb(panel.total_span)} total span, {len(panel.gene_names)} distinct labels.",
        f"Criteria bundle: `{report.bundle_id}` - {len(report.reaches)} criteria, "
        f"{len(report.small_variant_criteria)} of which need small-variant calling.",
        "",
        "This compares a criteria table against a panel design. It validates neither. The",
        "criteria are still an unverified model draft; the panel is still marked",
        "`AS_FUSION_PANEL_V1_UNCONFIRMED`. Research use only.",
        "",
        "## Small-variant criteria the panel cannot support",
        "",
    ]

    blocked = report.blocked
    if not blocked:
        lines += ["None. Every small-variant criterion names only genes the panel targets.", ""]
    else:
        lines += ["| Criterion | Status | In panel | Missing |", "|---|---|---|---|"]
        for item in blocked:
            lines.append(
                f"| {item.display_name} | `{item.status}` | "
                f"{', '.join(item.genes_in_panel) or '-'} | "
                f"{', '.join(item.genes_missing) or '-'} |"
            )
        lines.append("")

    reachable = [item for item in report.small_variant_criteria if item.reachable]
    lines += ["## Small-variant criteria the panel already supports", ""]
    if not reachable:
        lines += ["None.", ""]
    else:
        lines += ["| Criterion | Genes | Panel span |", "|---|---|---|"]
        for item in reachable:
            spans = "; ".join(
                f"{gene} {_kb(sum(t.span for t in panel.targets_for(gene)))}"
                for gene in item.genes_in_panel
            )
            lines.append(f"| {item.display_name} | {', '.join(item.genes_in_panel)} | {spans} |")
        lines.append("")
        lines += [
            "Span is the enriched territory carried in the BED, which includes the design's",
            "~10 kb flanks. It is not a claim that the gene's coding exons are covered; no",
            "coordinate here is checked against a gene model.",
            "",
        ]

    lines += [
        "## Genes a small-variant criterion needs and the panel does not target",
        "",
    ]
    genes = report.genes_to_add
    if not genes:
        lines += ["None.", ""]
    else:
        lines += [
            f"{len(genes)} genes: {', '.join(genes)}.",
            "",
            "Extending the design is a laboratory decision, not a software one. Adaptive",
            "sampling adds no reagent cost per region; the cost is that a larger enriched",
            "territory divides the same yield, so on-target depth falls.",
            "",
        ]

    return "\n".join(lines)


def _repo_root() -> Path:
    """Locate the checkout that carries the shipped panel and criteria.

    An editable install leaves ``__file__`` inside the checkout, which is how CI runs. A
    wheel does not, so fall back to the working directory, which is the convention the other
    tests here already follow.
    """
    candidate = Path(__file__).resolve().parents[2]
    if (candidate / "configs" / "panels").is_dir():
        return candidate
    return Path.cwd()


REPO_ROOT = _repo_root()
DEFAULT_PANEL = REPO_ROOT / "configs/panels/aml_fusion_adaptive_sampling.grch38.buffered.bed"
DEFAULT_CRITERIA = (
    REPO_ROOT / "configs/knowledge_bundles/GUIDELINE_CRITERIA_DRAFT_v0/guideline_criteria.v0.1.json"
)
GENERATED_REPORT = REPO_ROOT / "docs/PANEL_REACHABILITY.md"


def build_default_report() -> tuple[ReachabilityReport, Panel]:
    """Assess the shipped panel against the shipped criteria draft.

    Uses ``load_for_review`` because the whole point is to audit the draft *before* review;
    the reportable loader would refuse it, correctly.
    """
    from ontseq_platform.guideline_criteria import load_for_review

    bundle = load_for_review(DEFAULT_CRITERIA)
    panel = load_panel(DEFAULT_PANEL)
    return assess_panel(bundle.criteria, panel, bundle_id=bundle.bundle_id), panel


def main() -> int:
    report, panel = build_default_report()
    GENERATED_REPORT.write_text(format_report(report, panel) + "\n", encoding="utf-8")
    print(f"wrote {GENERATED_REPORT.relative_to(REPO_ROOT)}")
    print(
        f"{len(report.blocked)} of {len(report.small_variant_criteria)} small-variant criteria "
        f"are not supported by this panel; {len(report.genes_to_add)} genes are missing"
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
