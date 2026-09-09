ONTSeq Desktop v0.7.1 engineering shell

Research use only. Not clinically validated.
The executable starts and controls the existing ONTSeq backend; it does not replace analytical validation.
Active analysis profiles: AML_LCWGS_GRCh37, AML_AS_111_GRCh37,
AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25, AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25,
AML_LCWGS_GRCh38, AML_AS_111_GRCh38, AML_LCWGS_GRCh38_CANONICAL25 and
AML_AS_111_GRCh38_CANONICAL25. The unsuffixed profiles require their complete pinned assembly
dictionary; CANONICAL25 requires exactly chr1-22, chrX, chrY, chrM with build-specific lengths.
The GRCh37 Adaptive-Sampling profiles use only AML_AS_111_GRCh37_v1, never the GRCh38 panel:
110/111 selection intervals are mapped, 109 are reciprocal-exact and 108 have native GENCODE-19
ROIs; ACACA requires roundtrip review and CT45A2, IGH and GPR128 remain mapping/ROI gaps.
Default WSL resource root: ~/.local/share/ontseq/resources-v0.7.1 (under $HOME).
The default follows the software release; saved custom resource paths stay selected. GRCh37 and GRCh38
resources remain isolated by profile and are never mixed or lifted over automatically.
