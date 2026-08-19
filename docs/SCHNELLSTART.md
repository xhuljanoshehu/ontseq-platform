# Schnellstart — ONTSeq Platform

Diese Anleitung richtet sich an die Personen, die die Pipeline im Labor bedienen. Sie ist
bewusst auf Deutsch, während die übrige Dokumentation in `docs/` den technischen Nachweis
auf Englisch führt: das hier ist eine Bedienungsanleitung, kein Entwurfsdokument.

Sie ist auf eine konkrete Umgebung zugeschnitten: **GridION oder PromethION, Basecalling
läuft auf dem Gerät, Adaptive Sampling, Referenz `hg38` in UCSC-Benennung.**

---

## 0. Was Sie bekommen — und was nicht

**Nichts hiervon ist klinisch validiert.** Jeder Bericht trägt den Hinweis, dass er nur zu
Forschungszwecken dient, und jede Freigabe verlangt eine ausdrückliche Vidierung.

Ein Lauf liefert:

| Ergebnis | Zustand |
|---|---|
| Eingangsprüfung des BAM gegen die gelockte Referenz | läuft |
| QC-Kennzahlen (cramino) | läuft |
| Strukturvarianten (Sniffles2) | läuft |
| ISCN-Vorschlag aus den gefundenen Ereignissen | läuft |
| HTML-Bericht, Excel-Mappe, Ergebnis-JSON | läuft |
| Release-Bundle mit Prüfsummen | läuft |
| Vidierung mit Audit-Trail | läuft |
| **Kopienzahl (CNV)** | **`NOT_RUN` — kein Caller verdrahtet** |
| **Basecalling aus POD5** | **`NOT_RUN` — brauchen Sie nicht, das Gerät hat es getan** |
| **Wissensbasis-Annotation** | separater Befehl, ClinVar, siehe Schritt 9 |

`NOT_RUN` heißt: **diese Stufe wurde nicht ausgeführt**. Es ist kein negativer Befund. Der
Bericht sagt das an jeder Stelle ausdrücklich.

---

## 1. Voraussetzungen

- Linux, x86-64
- `micromamba` oder `conda` (holt die vier externen Werkzeuge aus bioconda — dafür wird
  einmalig Internet gebraucht)
- Plattenplatz: für die Auswertung eines Laufs mindestens **50 GB frei**, zusätzlich zu den
  Rohdaten. Ein vollständiger GridION-Lauf belegt inklusive POD5 rund **420 GB**.
- Die Rohdaten selbst müssen **nicht** kopiert werden; die Pipeline liest sie an Ort und
  Stelle und schreibt nur in ihr eigenes Ausgabeverzeichnis.

---

## 2. Auspacken und Umgebung anlegen

```bash
tar xzf ontseq-platform-<commit>.tar.gz
cd ontseq-platform

micromamba create -y -n ontseq -f workflow/envs/aligned_bam.yaml
micromamba activate ontseq
pip install -e .
```

Die Umgebung bringt `samtools 1.24`, `minimap2 2.28`, `sniffles 2.8.0` und `cramino 1.3.0`
in exakt den Versionen mit, gegen die die CI läuft. Sie liegt getrennt von allem, was sonst
auf der Maschine installiert ist — Ihre bestehende Illumina-Kette bleibt unberührt.

Prüfen, dass alles da ist:

```bash
ontseq --help
samtools --version | head -1
sniffles --version
```

---

## 2b. Für eine Maschine ohne Netz: der Container

Überspringen Sie diesen Abschnitt, wenn die Auswertungsmaschine ins Internet kann. Die
Pipeline selbst geht **nie** ins Netz; Schritt 2 braucht es einmalig, danach nie wieder.

Ist die Maschine dauerhaft getrennt, ist ein Container der saubere Weg. Alles steckt in
**einer Datei** — Werkzeuge, Bibliotheken, Betriebssystemebene —, und der Digest des
Abbilds nagelt alles zusammen fest.

### Auf einer Maschine **mit** Netz

```bash
docker build -f containers/Dockerfile -t ontseq:0.3.0 .
docker save ontseq:0.3.0 | gzip > ontseq-0.3.0.tar.gz
sha256sum ontseq-0.3.0.tar.gz > ontseq-0.3.0.tar.gz.sha256
```

Größenordnung: rund 0,5–1,5 GB. Bevorzugt Ihre Umgebung Apptainer — in klinischen
Rechenumgebungen die Regel —, liest das dasselbe Archiv:

```bash
apptainer build ontseq-0.3.0.sif docker-archive://ontseq-0.3.0.tar
```

### Der Transport

Ein **eigener Datenträger nur für diesen Zweck**, der nie etwas anderes transportiert und
nie ein Gerät außerhalb des Hauses berührt. Drei Regeln machen daraus eine
nachvollziehbare Übergabe statt eines Dateikopiervorgangs:

1. **Die Prüfsumme reist getrennt.** Nicht auf demselben Träger — auf Papier, per Mail, im
   LIMS. Eine Prüfsumme, die neben der Datei liegt, prüft nichts; sie bestätigt nur, dass
   die Datei zu sich selbst passt.
2. **Nur in eine Richtung.** Software hinein, nichts heraus. Ergebnisse verlassen die
   Maschine auf dem Weg, den Ihre Datenschutzregelung dafür vorsieht — nicht auf diesem
   Träger. Auf der Zielmaschine schreibgeschützt einhängen (`mount -o ro`), dann erzwingt
   das die Maschine statt der Disziplin.
3. **Beschriften.** Ein unbeschrifteter Datenträger in einer Schublade ist binnen eines
   Monats ein Allzweck-Datenträger, und die ganze Maßnahme ist dahin.

### Auf der Maschine **ohne** Netz

```bash
sha256sum -c ontseq-0.3.0.tar.gz.sha256      # zuerst, immer
gunzip -c ontseq-0.3.0.tar.gz | docker load

docker run --rm -v /daten:/daten -v /auswertung:/auswertung \
  ontseq:0.3.0 preflight /daten/probe.manifest.yaml \
  --reference-lock /daten/hg38.reference-lock.json --run-id RUN_001
```

Die Datenverzeichnisse werden hineingereicht und mit **absoluten** Pfaden angesprochen. Alle
weiteren Schritte dieser Anleitung funktionieren unverändert, nur mit `docker run ...` statt
`ontseq ...` davor.

**Heben Sie den Datenträger auf.** Er trägt dann genau die Fassung, die installiert wurde,
mit ihrer Prüfsumme — und ist damit Teil des Nachweises, nicht nur ein Transportmittel.

> Das Container-Rezept wird bei jedem CI-Lauf gebaut, die vier Werkzeugversionen darin
> werden geprüft, und die Pipeline läuft im Abbild einmal komplett durch. Es ist also
> nachgewiesen, nicht behauptet.

---

## 2c. Die Oberfläche im Browser

Statt der Kommandozeile geht es auch per Browser. Der Dienst startet die Läufe, die
Oberfläche zeigt sie an.

```bash
ontseq serve \
  --reference-lock /pfad/zu/hg38.reference-lock.json \
  --allow-root /mnt/p/NANOPORE \
  --output-dir /mnt/p/AUSWERTUNG
```

Danach öffnet sich `http://127.0.0.1:8765/`. Unter Windows mit WSL2 genügt eine
`ONTSeq.bat` auf dem Desktop:

```bat
@echo off
start "" http://127.0.0.1:8765/
wsl -d Ubuntu -- bash -lc "ontseq serve --reference-lock ... --allow-root /mnt/p/NANOPORE"
```

Doppelklick, Browser geht auf, Oberfläche da.

### Was der Dienst nicht tut

Er rechnet nichts. Er ruft dieselbe Funktion auf, die auch `ontseq run` aufruft, liest den
Fortschritt aus derselben `provenance/run.json`, die auch `ontseq status` liest, und reicht
die Ergebnisse durch. Ein Fehler dort kann ein Ergebnis verlieren — verändern kann er keines.

### Sicherheit, in vier Sätzen

- **Nur `127.0.0.1`.** Kein Netz, keine anderen Rechner.
- **Token pro Prozess.** Die Seite bekommt es beim Ausliefern eingesetzt. Eine von der
  Platte geöffnete Kopie der Datei hat keins und kommt deshalb nicht durch — beabsichtigt.
- **`--allow-root` ist eine Whitelist.** Was nicht darunter liegt, ist nicht lesbar, egal
  welcher Pfad angefragt wird. Symlinks werden vorher aufgelöst.
- **`Host` und `Origin` werden geprüft.** Das ist der Schutz gegen DNS-Rebinding und gegen
  eine andere offene Registerkarte, die im Hintergrund Läufe startet.

Wird `--allow-root` weggelassen, ist der Dienst nicht etwa offen, sondern liest gar nichts.

### Der Dateidialog kommt vom Dienst, nicht vom Browser

Ein Browser gibt JavaScript nur den Dateinamen, nie den Pfad. Der Dienst listet deshalb die
freigegebenen Verzeichnisse selbst auf. Nebeneffekt: die Liste zeigt gleich, ob zu einer BAM
ein `.bai` existiert — statt dass der Lauf Stunden später daran scheitert.

### Vier Zustände, nicht zwei

Die Fortschrittsanzeige unterscheidet **✓ abgeschlossen**, **◇ kein Nachweis**,
**✕ fehlgeschlagen** und **– nicht ausgeführt**. Das ist kein Detail: „nicht ausgeführt" und
„nichts gefunden" bedeuten das Gegenteil voneinander, und CNV steht heute auf *nicht
ausgeführt*. Eine Oberfläche, die beides als leeren Kreis malt, macht genau den Fehler, den
der Rest des Systems verhindert.

**Vidiert wird nicht in der Oberfläche.** Die Freigabe läuft über `ontseq review record`.

---

## 3. Referenz einmalig locken

Das machen Sie **einmal pro Referenz**, nicht pro Lauf. Der Lock hält fest, welche Contigs
die Referenz hat und wie ihr Index aussieht; ab dann prüft die Pipeline bei jedem Lauf, ob
sich daran etwas geändert hat.

```bash
ontseq reference-lock \
  --fai /pfad/zu/hg38.fa.fai \
  --reference-id UCSC_HG38_CHR_ONLY_2024-11 \
  --genome-build GRCh38 \
  --output /pfad/zu/hg38.reference-lock.json
```

**Zur `--reference-id`:** wählen Sie einen Namen, der die Herkunft trägt — Quelle und
Stand. `hg38` allein genügt nicht; es gibt zu viele Dateien, die so heißen.

> **Hinweis zur PAR.** Ist in Ihrer Referenz die pseudoautosomale Region auf chrY nicht
> hart maskiert, verlieren Reads in PAR1 auf chrX ihre Mapping-Qualität und fallen aus der
> Abdeckung. Falls Sie dort Ziele haben, prüfen Sie das einmal — es betrifft dann jeden
> Lauf gleichermaßen. Siehe `docs/PIPELINE_EXECUTION.md` §9.

---

## 4. Manifest für eine Probe schreiben

Das Manifest ist die einzige Stelle, an der Sie sagen, **was** analysiert wird. Vorlage:
`examples/manifests/gridion_adaptive_sampling.example.yaml`. Kopieren und anpassen:

```yaml
schema_version: 0.1.0
sample_id: PROBE_001            # pseudonymisiert, keine Klarnamen
run_id: 260611_RAD114_AS_PROBE  # frei wählbar, taucht im Ausgabepfad auf

input:
  kind: aligned_bam
  path: /daten/260611_.../sample.bam
  index_path: /daten/260611_.../sample.bam.bai

assay:
  mode: adaptive_sampling
  genome_build: GRCh38
  reference_id: UCSC_HG38_CHR_ONLY_2024-11    # identisch zum Lock
  target_bed: /data/adaptive_sampling/250611_fusion_panel_with_buffer.bed
  target_bed_version: "250611"

analysis:
  profile: adaptive_sampling
  modules: [qc, sv, iscn, report]
  intent: somatic               # AML fragt nach erworbenen Veränderungen

privacy:
  pseudonymized: true
  contains_direct_identifiers: false
  cloud_upload_approved: false
```

Drei Dinge, die die Pipeline erzwingt:

- **`reference_id` muss zum Lock passen.** Sonst bricht die Eingangsprüfung ab. Das ist der
  Schutz davor, ein BAM gegen die falsche Referenz auszuwerten.
- **`adaptive_sampling` verlangt `target_bed` und `target_bed_version`.** Ohne Angabe, was
  angereichert wurde, ist eine Abdeckungszahl nicht interpretierbar.
- **`contains_direct_identifiers: true` wird abgelehnt.** Keine Klarnamen, keine
  Geburtsdaten, keine Fallnummern im Manifest.

`intent: somatic` hat keinen Vorgabewert und muss gesetzt werden, sonst meldet jede spätere
ClinVar-Annotation ihren Geltungsbereich als `unknown` — richtig, aber nutzlos.

---

## 5. Preflight — prüfen, bevor Sie starten

```bash
ontseq preflight probe_001.manifest.yaml \
  --reference-lock /pfad/zu/hg38.reference-lock.json \
  --run-id 260611_RAD114_AS_PROBE \
  --require-free-gb 50
```

Prüft in Sekunden und **ohne jede Nebenwirkung** — es wird kein Verzeichnis angelegt, keine
Sperre gesetzt, nichts geschrieben: Existiert das BAM? Passt der Index? Stimmen die Contigs
mit dem Lock überein? Sind alle Werkzeuge da und in der richtigen Version? Ist genug Platz
frei?

| Anzeige | Bedeutung |
|---|---|
| `ok` | geprüft, Voraussetzung erfüllt |
| `FAIL` | geprüft, der Lauf kann so nicht gelingen — **einziger blockierender Zustand** |
| `warn` | der Lauf kann laufen, jemand sollte es trotzdem wissen |
| `????` | von hier aus nicht feststellbar |
| `--` | trifft auf diese Eingabeart nicht zu |

Exit 0 = nichts blockiert, Exit 2 = mindestens eine Voraussetzung fehlt.

---

## 6. Den Lauf starten

```bash
ontseq run probe_001.manifest.yaml \
  --reference-lock /pfad/zu/hg38.reference-lock.json \
  --run-id 260611_RAD114_AS_PROBE \
  --output-dir /auswertung/runs \
  --threads 8
```

Ergebnis landet unter `/auswertung/runs/260611_RAD114_AS_PROBE/PROBE_001/`.

**Abgebrochen? Einfach denselben Befehl noch einmal.** Fertige Stufen werden nicht
wiederholt — die Pipeline prüft Prüfsummen, Parameter und Werkzeugversionen und setzt genau
dort wieder an, wo sie aufgehört hat. Sie erkennt auch, wenn sich zwischendurch etwas
geändert hat, und rechnet das Betroffene dann neu.

Ein Envelope kann immer nur von **einem** Lauf gleichzeitig bearbeitet werden. Ein zweiter
Aufruf bricht mit Exit 4 ab und nennt, wer die Sperre hält.

---

## 7. Ergebnis ansehen

```bash
ontseq status --output-dir /auswertung/runs
ontseq status --output-dir /auswertung/runs --verbose      # je Stufe
```

Im Envelope liegen:

| Pfad | Inhalt |
|---|---|
| `reports/` | HTML-Bericht und Excel-Mappe |
| `normalized/` | Ergebnis als JSON (maschinenlesbar) |
| `qc/`, `evidence/` | QC-Kennzahlen, VCF |
| `provenance/run.json` | jede Stufe, jede Werkzeugversion, jeder Parameter |
| `release/` | die freizugebenden Dateien plus `checksums.sha256` |

Prüfsummen unabhängig nachrechnen:

```bash
cd /auswertung/runs/260611_RAD114_AS_PROBE/PROBE_001/release
sha256sum -c checksums.sha256
```

---

## 8. Vidieren

```bash
# Zustand ansehen
ontseq review status /auswertung/runs/260611_RAD114_AS_PROBE/PROBE_001

# freigeben
ontseq review record /auswertung/runs/260611_RAD114_AS_PROBE/PROBE_001 \
  --decision accepted \
  --reviewer "Dr. Muster" \
  --note "Befund geprüft, ISCN-Vorschlag übernommen"

# oder ablehnen
ontseq review record ... --decision rejected --reviewer "..." --note "Grund"
```

Zwei Eigenschaften, die Sie kennen sollten:

- **Die Vidierung bindet an den Inhalt, nicht an den Pfad.** Ändert sich die freigegebene
  Datei nachträglich, wird die Freigabe als `STALE` gemeldet — sie gilt dann nicht mehr.
- **Ein vidierter Envelope kann nicht erneut durchlaufen werden.** Der Versuch bricht mit
  Exit 7 ab. Das lässt sich nicht per Flag übergehen: es soll nicht möglich sein, das
  nachträglich zu ändern, was jemand unterschrieben hat.

Der Prüfpfad ist eine Hash-Kette. Eine nachträgliche Änderung ist damit **erkennbar** —
aber nicht **verhinderbar**: es gibt keinen Schlüssel. Wer die Datei schreiben darf, kann
die ganze Kette konsistent neu schreiben. Für eine echte Signatur fehlt die
Schlüsselinfrastruktur.

---

## 9. Optional: ClinVar-Annotation

```bash
ontseq annotate /auswertung/runs/.../normalized/PROBE_001.result.json \
  --clinvar /pfad/zu/variant_summary.txt \
  --release 2026-08-01 \
  --output PROBE_001.annotated.json
```

**Lesen Sie den Hinweis im Bericht.** ClinVar klassifiziert nach ACMG-Keimbahnregeln; eine
AML-Abklärung stellt eine somatische Frage. Ein Eintrag, dessen Herkunft nicht zur
Fragestellung passt, wird **behalten und markiert**, nie herausgefiltert — und keine
Annotation macht irgendetwas berichtsfähig. Im Excel-Blatt `11_Annotations` steht die
Leseregel in Zeile 1, und nicht passende Zeilen sind farbig hinterlegt.

---

## 10. Optional: Watch-Ordner (unbeaufsichtigt)

```bash
ontseq watch --config watch.yaml
```

Der Watcher nimmt einen Lauf erst auf, wenn MinKNOW die Datei `final_summary_*.txt`
geschrieben hat, und führt jede Probe **höchstens einmal** aus. Fehlgeschlagene Proben
werden **nicht** automatisch wiederholt — dafür gibt es `--retry-failed`, ausdrücklich.

Der Watcher rät nichts. Was er nicht aus dem Ordner ableiten kann, muss in der
Konfiguration stehen. Details: `docs/PIPELINE_EXECUTION.md` §7.

---

## 11. Wenn etwas schiefgeht

| Exit-Code | Bedeutung |
|---|---|
| 0 | alles in Ordnung |
| 2 | Lauf oder Probe fehlgeschlagen, oder eine Voraussetzung fehlt |
| 3 | Karyotyp nur teilweise umgewandelt |
| 4 | Envelope ist gesperrt — ein anderer Lauf arbeitet daran |
| 5 | Watcher-Konfiguration unbrauchbar |
| 6 | abgebrochen oder unvollendet |
| 7 | Envelope trägt eine gültige Vidierung |

Ein abgebrochener Lauf ist **kein** Datenverlust: jeder Schreibvorgang ist atomar, es liegt
immer entweder die vorige vollständige Datei oder gar keine. Erneut starten genügt.

---

## 12. Was bewusst fehlt

Damit niemand danach sucht:

- **Kein CNV-Caller.** Die Stufe ist im Graphen deklariert und meldet `NOT_RUN`. Welcher
  Caller dort eingehängt wird, ist eine fachliche Entscheidung.
- **Keine somatische Wissensquelle.** ClinVar ist eine Keimbahn-Quelle. OncoKB, CIViC,
  COSMIC oder ELN/ICC-Kriterien sind nicht eingebunden.
- **Keine Signatur** am Release-Bundle. `signature_status` steht wörtlich auf `unsigned`.
- **Keine Interpretation.** Die Pipeline endet beim geprüften Befund. Die Bewertung ist
  ärztliche Aufgabe, und dafür ist der Bericht gemacht.
