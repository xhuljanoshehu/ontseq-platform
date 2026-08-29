ONTSeq Desktop v0.5.3 engineering shell

Research use only. Not clinically validated.
The executable starts and controls the existing ONTSeq backend; it does not replace analytical validation.
Active analysis profiles: AML_LCWGS_GRCh38, AML_AS_111_GRCh38,
AML_LCWGS_GRCh38_CANONICAL25 and AML_AS_111_GRCh38_CANONICAL25.
The unsuffixed profiles require the complete Primary-Assembly dictionary; CANONICAL25 requires
exactly chr1-22, chrX, chrY, chrM. All four use the same installed GRCh38 bundles.
Default WSL resource root: ~/.local/share/ontseq/resources (under $HOME). GRCh37 resources are not mixed into these profiles.
