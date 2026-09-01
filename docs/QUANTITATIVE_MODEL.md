# Quantitative model: tumour fraction, copy-number baseline and subclones

GENERATED FILE - do not edit by hand. Regenerate with:

    python -m ontseq_platform.quantitation

Every number below is computed by `ontseq_platform.quantitation` at generation time,
so the document cannot drift from the code. Nothing here is validated for clinical
use, and no value is an assay adequacy threshold. Research use only.

## The problem

A sample is a mixture of tumour and normal cells. A copy-number ratio, a methylation
beta value and an allele fraction are all diluted by the normal fraction, and none of
them states by how much. A copy-number call reported without the tumour fraction is a
ratio presented as a count.

## Why the estimator cannot come from copy number

The tumour fraction is what converts an observed ratio into an integer copy number, so
deriving it from those same ratios is circular. The anchor has to be independent of
copy number: the allele fraction of a clonal, heterozygous somatic SNV in a
copy-neutral diploid region. Half the alleles in a tumour cell carry it and none in a
normal cell do, so `VAF = f / 2`, and therefore `f = 2 x VAF`.

The identity holds *only* where the region is genuinely copy-neutral and diploid. Under
a deletion, an amplification or copy-neutral LOH the allele fraction shifts and the
resulting tumour fraction is wrong invisibly. Copy-number calls select the regions the
estimate may be taken from; they never supply the number. `copy_number_state` is
therefore a required argument with no default.

## What follows once the tumour fraction is known

1. **Copy-number baseline.** An observed ratio resolves to a tumour copy number.
2. **Methylation.** Beta values can be corrected for normal-cell dilution.
3. **Subclones.** A second variant's allele fraction becomes a cancer cell fraction,
   `CCF = 2 x VAF / f` - the proportion of the tumour carrying it.

## Why a fixed copy-ratio threshold is unsafe

The same true copy number produces a different observed ratio at every tumour fraction.
A threshold chosen for pure tumour misses real events in a diluted sample.

| True copy number | f = 20% | f = 30% | f = 50% | f = 70% | f = 90% |
|---|---|---|---|---|---|
| 0 | 0.80 | 0.70 | 0.50 | 0.30 | 0.10 |
| 1 | 0.90 | 0.85 | 0.75 | 0.65 | 0.55 |
| 2 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3 | 1.10 | 1.15 | 1.25 | 1.35 | 1.45 |
| 4 | 1.20 | 1.30 | 1.50 | 1.70 | 1.90 |

A one-copy deletion shows a ratio of 0.75 at half tumour content, not 0.50. A filter
set at a fixed distance from 1.0 is a filter whose meaning changes with every sample.

## Detection limits by depth

Per-base error rate 1% (technical default), alpha 0.01,
power 0.95. `Min reads` is the fewest variant reads that sequencing error
alone would produce with probability at most alpha - the floor any real variant must
clear. It is a property of depth and error rate, not of a caller: a caller with a
better error model does better than this, never worse.

| Depth | Min variant reads | Smallest detectable VAF |
|---|---|---|
| 9x (measured off-target) | 2 | 42.9% |
| 20x | 3 | 28.3% |
| 30x | 3 | 19.5% |
| 50x | 4 | 14.8% |
| 80x (measured on-target) | 4 | 9.4% |
| 100x | 5 | 8.9% |
| 200x | 7 | 5.8% |
| 500x | 12 | 3.6% |

## Smallest resolvable subclone

The question a haematologist asks. Two limits stack: the allele fraction has to clear
the error floor, and the tumour has to be a large enough part of the sample for that
allele fraction to correspond to a small subclone. Values are cancer cell fractions -
the proportion of *tumour* cells carrying the variant.

| Depth | f = 20% | f = 30% | f = 50% | f = 70% | f = 90% |
|---|---|---|---|---|---|
| 9x | not resolvable | not resolvable | not resolvable | not resolvable | 95.4% |
| 20x | not resolvable | not resolvable | not resolvable | 80.7% | 62.8% |
| 30x | not resolvable | not resolvable | 78.1% | 55.8% | 43.4% |
| 50x | not resolvable | 98.6% | 59.1% | 42.2% | 32.9% |
| 80x | 94.1% | 62.7% | 37.6% | 26.9% | 20.9% |
| 100x | 89.2% | 59.5% | 35.7% | 25.5% | 19.8% |
| 200x | 58.3% | 38.9% | 23.3% | 16.7% | 13.0% |
| 500x | 36.2% | 24.1% | 14.5% | 10.3% | 8.0% |

Read against this assay's measured depths:

* **On-target (80x)** resolves allele fractions down to about 9.4%. At a tumour fraction of 50% that is a subclone of about 37.6% of tumour cells - major subclones only, not minor ones.
* **Off-target (9x)** resolves nothing below 95.4%, and only at the highest tumour
  fraction in this table; below that it resolves no subclone at all. That is a
  property of the depth, not of the caller, and no choice of software changes it.

So the entire quantitative model lives inside the panel. Genome-wide, this assay can
carry copy number and methylation, and cannot carry allele-fraction quantitation.

## What is not implemented

**All of the above consumes allele counts, and ONTSeq has no small-variant caller
wired in.** There is no stage, no pinned caller, and no validation. This module is the
arithmetic and its limits, written now because both are decidable now; nothing in the
pipeline calls it. The dependency is the same one that blocks seven of the
twenty-four drafted guideline criteria.

Also not addressed here, and each a separate problem:

* Distinguishing somatic from germline without a matched normal. A germline
  heterozygous SNP sits at VAF 0.5 whatever the tumour content and would report a
  tumour fraction of 1. The module refuses that value rather than returning it, which
  is a guard, not a solution.
* Clonality. The module cannot tell a clonal variant from a subclonal one without
  being told the tumour fraction, and the tumour fraction is what a clonal variant is
  used to establish. Breaking that requires either a known-clonal marker or joint
  estimation over many sites.
* Multi-site inference. Real purity and subclone calling fits a mixture model over
  many variants at once. These functions are per-site and deliberately so.

## In-silico dilution

`expected_vaf` is the forward model, and it is the generator for a dilution series:
choose a tumour fraction and a subclone size, and it gives the allele fraction the
mixed reads should carry. Mixing reads from a known-positive sample with reads from a
normal at defined ratios measures the sensitivity curve against a known truth, needs
no new patient material and no new sequencing run, and validates the table above
rather than trusting it. The table is a ceiling derived from counting statistics; a
real caller will do worse, and the gap between them is the number worth having.

