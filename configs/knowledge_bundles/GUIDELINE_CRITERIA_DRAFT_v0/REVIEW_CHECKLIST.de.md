# Prüfliste: Leitlinienkriterien ELN 2022 / WHO 2022 / ICC 2022

> **Dieser Entwurf ist ungeprüft.** Jeder Eintrag wurde von einem Sprachmodell aus dem
> Gedächtnis geschrieben und *nicht* gegen die Leitlinientexte abgeglichen. Einträge
> können falsch, unvollständig, veraltet oder erfunden sein.

## So wird geprüft

1. Leitlinientext neben diese Liste legen.
2. Je Eintrag: Wortlaut vergleichen, `guideline_reference` mit Abschnitt/Tabelle füllen.
3. `verification` von `unverified_model_draft` auf `verified` oder `rejected` setzen.
4. Am Ende `provenance.reviewer` und `provenance.review_date` eintragen.

Solange ein Eintrag `unverified_model_draft` trägt, darf ihn keine Software verwenden.

**24 Einträge**, davon 14 mit dem
heutigen Assay überhaupt bestimmbar.

## Zuerst prüfen — hier vermute ich Änderungen gegenüber ELN 2017

- [ ] **Mutated NPM1 without FLT3-ITD** (`ELN2022-FAV-NPM1-NO-FLT3ITD`)
      Vermutete Änderung: allelic-ratio handling believed changed in 2022
      VERIFY CAREFULLY: ELN 2017 used the FLT3-ITD allelic ratio; I believe 2022 dropped the ratio. Confirm the exact 2022 condition.

- [ ] **In-frame bZIP mutated CEBPA** (`ELN2022-FAV-CEBPA-BZIP`)
      Vermutete Änderung: biallelic -> in-frame bZIP believed changed in 2022
      VERIFY CAREFULLY: I believe 2022 changed this from biallelic CEBPA to in-frame bZIP. Confirm the exact criterion.

- [ ] **t(8;16)(p11;p13) KAT6A::CREBBP** (`ELN2022-ADV-KAT6A-CREBBP`)
      Vermutete Änderung: believed added in 2022
      I believe this was ADDED in 2022. Confirm presence and bands.

- [ ] **t(3q26.2;v) MECOM(EVI1)-rearranged** (`ELN2022-ADV-MECOM-REARRANGED`)
      Vermutete Änderung: believed added in 2022
      I believe the generalized MECOM-rearranged criterion was ADDED in 2022.

- [ ] **Mutated ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1 and/or ZRSR2** (`ELN2022-ADV-MDS-RELATED-GENES`)
      Vermutete Änderung: believed added in 2022
      VERIFY THE GENE LIST MEMBER BY MEMBER. I believe this set was added in 2022 and that it does not apply to otherwise favourable-risk AML. Confirm both.

## Vollständige Liste

### ELN 2022 — günstig

- [ ] **t(8;21)(q22;q22.1) RUNX1::RUNX1T1** — berechenbar
      `ELN2022-FAV-RUNX1-RUNX1T1` · ISCN `t(8;21)(q22;q22.1)`
      Confirm the ELN 2022 wording and whether any blast or qualifier condition applies.

- [ ] **inv(16)(p13.1q22) / t(16;16)(p13.1;q22) CBFB::MYH11** — berechenbar
      `ELN2022-FAV-CBFB-MYH11` · ISCN `inv(16)(p13.1q22)`
      Confirm both the inv(16) and t(16;16) forms are listed together.

- [ ] **Mutated NPM1 without FLT3-ITD** — **braucht Variantencalling — fehlt**
      `ELN2022-FAV-NPM1-NO-FLT3ITD`
      VERIFY CAREFULLY: ELN 2017 used the FLT3-ITD allelic ratio; I believe 2022 dropped the ratio. Confirm the exact 2022 condition.

- [ ] **In-frame bZIP mutated CEBPA** — **braucht Variantencalling — fehlt**
      `ELN2022-FAV-CEBPA-BZIP`
      VERIFY CAREFULLY: I believe 2022 changed this from biallelic CEBPA to in-frame bZIP. Confirm the exact criterion.

### ELN 2022 — intermediär

- [ ] **Mutated NPM1 with FLT3-ITD** — **braucht Variantencalling — fehlt**
      `ELN2022-INT-NPM1-WITH-FLT3ITD`
      Confirm placement in intermediate regardless of allelic ratio.

- [ ] **Wild-type NPM1 with FLT3-ITD, without adverse-risk lesions** — **braucht Variantencalling — fehlt**
      `ELN2022-INT-FLT3ITD-WT-NPM1`
      Confirm the exclusion clause for co-occurring adverse lesions.

- [ ] **t(9;11)(p21.3;q23.3) MLLT3::KMT2A** — berechenbar
      `ELN2022-INT-MLLT3-KMT2A` · ISCN `t(9;11)(p21.3;q23.3)`
      Confirm precedence: I believe t(9;11) takes priority over co-occurring adverse lesions.

- [ ] **Abnormalities not classified as favourable or adverse** — **braucht Variantencalling — fehlt**
      `ELN2022-INT-UNCLASSIFIED`
      This is the fallback bucket. It can only be applied once every other criterion has been evaluated, which this assay cannot do without variant calling.

### ELN 2022 — ungünstig

- [ ] **t(6;9)(p23;q34.1) DEK::NUP214** — berechenbar
      `ELN2022-ADV-DEK-NUP214` · ISCN `t(6;9)(p23;q34.1)`
      Confirm band designations.

- [ ] **t(v;11q23.3) KMT2A-rearranged** — berechenbar
      `ELN2022-ADV-KMT2A-REARRANGED` · ISCN `t(v;11q23.3)`
      IMPORTANT: confirm that t(9;11) is excluded here and sits in intermediate instead.

- [ ] **t(9;22)(q34.1;q11.2) BCR::ABL1** — berechenbar
      `ELN2022-ADV-BCR-ABL1` · ISCN `t(9;22)(q34.1;q11.2)`
      Confirm inclusion in ELN 2022 adverse.

- [ ] **t(8;16)(p11;p13) KAT6A::CREBBP** — berechenbar
      `ELN2022-ADV-KAT6A-CREBBP` · ISCN `t(8;16)(p11;p13)`
      I believe this was ADDED in 2022. Confirm presence and bands.

- [ ] **inv(3)(q21.3q26.2) / t(3;3)(q21.3;q26.2) GATA2, MECOM(EVI1)** — berechenbar
      `ELN2022-ADV-MECOM-INV3` · ISCN `inv(3)(q21.3q26.2)`
      Confirm both forms and the GATA2/MECOM annotation.

- [ ] **t(3q26.2;v) MECOM(EVI1)-rearranged** — berechenbar
      `ELN2022-ADV-MECOM-REARRANGED` · ISCN `t(3q26.2;v)`
      I believe the generalized MECOM-rearranged criterion was ADDED in 2022.

- [ ] **-5 or del(5q)** — berechenbar
      `ELN2022-ADV-MINUS5-DEL5Q` · ISCN `-5 / del(5q)`
      Confirm whether a minimal deleted region or size threshold is specified.

- [ ] **-7** — berechenbar
      `ELN2022-ADV-MINUS7` · ISCN `-7`
      Confirm whether del(7q) is included or only whole -7.

- [ ] **-17 or abn(17p)** — berechenbar
      `ELN2022-ADV-MINUS17-ABN17P` · ISCN `-17 / abn(17p)`
      Confirm the definition of abn(17p) and its relation to TP53.

- [ ] **Complex karyotype (>=3 unrelated abnormalities)** — **nur Untergrenze — Assay unvollständig**
      `ELN2022-ADV-COMPLEX-KARYOTYPE` · ISCN `>=3 abnormalities`
      CRITICAL: confirm the exact definition and its exclusions. This assay cannot count balanced rearrangements outside the panel, so the count is a LOWER BOUND and must never be reported as a negative complex-karyotype result.

- [ ] **Monosomal karyotype** — **nur Untergrenze — Assay unvollständig**
      `ELN2022-ADV-MONOSOMAL-KARYOTYPE`
      Confirm the definition and whether ELN 2022 still lists it. Same lower-bound limitation as complex karyotype.

- [ ] **Mutated TP53** — **braucht Variantencalling — fehlt**
      `ELN2022-ADV-TP53`
      Confirm any VAF threshold and the relationship to 17p loss.

- [ ] **Mutated ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1 and/or ZRSR2** — **braucht Variantencalling — fehlt**
      `ELN2022-ADV-MDS-RELATED-GENES`
      VERIFY THE GENE LIST MEMBER BY MEMBER. I believe this set was added in 2022 and that it does not apply to otherwise favourable-risk AML. Confirm both.

### WHO 2022 / ICC 2022 — definierende Aberrationen

- [ ] **AML with NUP98 rearrangement** — berechenbar
      `WHOICC-2022-NUP98-REARRANGED`
      Confirm presence in WHO 5th ed and in ICC, and note where the two differ.

- [ ] **AML with RBM15::MRTFA** — berechenbar
      `WHOICC-2022-RBM15-MRTFA`
      Note the MRTFA/MKL1 synonym; confirm which symbol each classification uses.

- [ ] **Blast percentage requirements differ between WHO 2022 and ICC 2022** — **nicht sequenzierbar**
      `WHOICC-2022-BLAST-THRESHOLD`
      CRITICAL DIVERGENCE: WHO 5th ed and ICC 2022 handle blast thresholds differently for genetically defined entities. Blast count is not a sequencing observable at all. Record explicitly which classification the report follows.

## Was der Assay strukturell nicht kann

Diese Grenzen verschwinden nicht durch Prüfen der Kriterien:

- **Kleinvarianten werden nicht gerufen.** NPM1, FLT3-ITD, CEBPA, TP53 und die
  MDS-assoziierten Gene sind damit heute unbestimmbar — nicht negativ, sondern
  ungemessen.
- **Komplexer und monosomaler Karyotyp** können nur als *Untergrenze* gezählt werden.
  Balancierte Umbauten außerhalb des Panels sieht das Verfahren nicht.
- **Blastenanteil** ist keine Sequenzierbeobachtung.

Solange das gilt, darf aus den übrigen Kriterien **keine ELN-Risikogruppe** abgeleitet
werden. Die korrekte Ausgabe ist „nicht bestimmbar“ mit Angabe, welche Kriterien fehlen.
