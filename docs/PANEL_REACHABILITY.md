# Panel reachability of the guideline criteria

GENERATED FILE - do not edit by hand. Regenerate with:

    python -m ontseq_platform.panel_reachability

Panel: `aml_fusion_adaptive_sampling.grch38.buffered.bed` - 111 targets, 17.03 Mb total span, 111 distinct labels.
Criteria bundle: `GUIDELINE_CRITERIA_DRAFT_v0` - 24 criteria, 7 of which need small-variant calling.

This compares a criteria table against a panel design. It validates neither. The
criteria are still an unverified model draft; the panel is still marked
`AS_FUSION_PANEL_V1_UNCONFIRMED`. Research use only.

## Small-variant criteria the panel cannot support

| Criterion | Status | In panel | Missing |
|---|---|---|---|
| Mutated ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1 and/or ZRSR2 | `some_named_genes_missing_from_panel` | RUNX1 | ASXL1, BCOR, EZH2, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2 |
| Mutated TP53 | `no_named_gene_in_panel` | - | TP53 |
| In-frame bZIP mutated CEBPA | `no_named_gene_in_panel` | - | CEBPA |
| Abnormalities not classified as favourable or adverse | `criterion_names_no_genes` | - | - |

## Small-variant criteria the panel already supports

| Criterion | Genes | Panel span |
|---|---|---|
| Mutated NPM1 without FLT3-ITD | NPM1, FLT3 | NPM1 43.1 kb; FLT3 117.3 kb |
| Mutated NPM1 with FLT3-ITD | NPM1, FLT3 | NPM1 43.1 kb; FLT3 117.3 kb |
| Wild-type NPM1 with FLT3-ITD, without adverse-risk lesions | NPM1, FLT3 | NPM1 43.1 kb; FLT3 117.3 kb |

Span is the enriched territory carried in the BED, which includes the design's
~10 kb flanks. It is not a claim that the gene's coding exons are covered; no
coordinate here is checked against a gene model.

## Genes a small-variant criterion needs and the panel does not target

10 genes: ASXL1, BCOR, CEBPA, EZH2, SF3B1, SRSF2, STAG2, TP53, U2AF1, ZRSR2.

Extending the design is a laboratory decision, not a software one. Adaptive
sampling adds no reagent cost per region; the cost is that a larger enriched
territory divides the same yield, so on-target depth falls.

