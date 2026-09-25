# MARLIN im automatischen ONTSeq-Analyselauf

Status: vom Nutzer am 2026-09-23 zur Umsetzung freigegeben; technische Abnahme ausstehend.
Datum: 2026-09-23

## Ziel und Erfolgskriterien

Wenn der Nutzer im Desktop die Methylierung auswählt, führt ONTSeq auch MARLIN aus.
Der normale Analysestart genügt; es sind keine separaten Terminalbefehle erforderlich.
Vor dem Start zeigt die Anwendung, ob die MARLIN-Voraussetzungen erfüllt sind.
Der Befund enthält den tatsächlichen Ausführungsstatus, die führende Modellklasse,
den Score, die Anzahl beobachteter Modell-CpGs und nachvollziehbare technische Nachweise.
Ein fehlendes Ergebnis wird niemals durch Beispielwerte oder einen festen Platzhalter ersetzt.

Der Auftrag umfasst Backend, Desktop/Browser, HTML/JSON/XLSX, Tests, Dokumentation,
Review, Integration in main und die Aktualisierung des vorhandenen lokalen Startprogramms.
Eine fachliche oder klinische Validierung wird damit nicht behauptet.

## Tatsächlich geprüfter Ausgangszustand

- Quellstand: ca7949728893c36e3ba541281dc1cbaed239ddc4.
- Der normale Ablauf hat einen Methylierungsschritt, aber keinen MARLIN-Schritt.
- `web/src/befund/main.jsx` zeigt im MARLIN-Abschnitt immer denselben Fehltext.
- Bestehende `marlin_*`-Module sind separate Forschungsbefehle mit einem R-Backend.
- Ihre bisherigen v1-Verträge sind auf GRCh37 begrenzt. Der native Übergang verlangt
  einen `MarlinBridgeLock` mit einem echten Vergleich derselben Probe.
- Das aktive Desktopprofil ist `AML_AS_111_GRCh38_CANONICAL25`.
- Der aktuelle offizielle MARLIN-Quellstand
  `442aa603415a54f62e7367794f9a31c6bc20fc2d` veröffentlicht hg19- und hg38-Sondenkarten.
  Die hg19-Begrenzung ist daher eine Einschränkung unseres bestehenden Adapters.
- Auf diesem Rechner sind Originalmodell, Merkmalsreihenfolge, Klassenannotationen,
  hg19/hg38-Sondenkarten und eine separat paketierte TensorFlow-CPU-Laufzeit vorhanden.
- Frühere tatsächliche Modellaufrufe und Offline-Wiederherstellung wurden mit künstlichen
  Daten geprüft. Dies ist kein biologischer Vergleich derselben Probe und erfüllt den
  bestehenden `validated_same_specimen_bridge`-Nachweis nicht.

## Gewählter technischer Ansatz

Ein eigener, ausdrücklich als nicht analytisch validiert gekennzeichneter Forschungsadapter
verbindet den aktuellen Analyselauf mit der separat installierten MARLIN-Laufzeit.
Vorhandene strengere v1-Befehle und ihre Vergleichsnachweise bleiben in ihrer Bedeutung erhalten.
Der neue Adapter erhält eine eigene Versionskennung und eigene Ressourcen-/Laufzeitnachweise;
er erzeugt keinen fiktiven v1-Bridge-Lock.

Direkte Auswertung mit der zum BAM passenden offiziellen Sondenkarte:
GRCh38 nutzt hg38, GRCh37 nutzt hg19. Eine zusätzliche Ausrichtung oder ein automatisches
Liftover ist nicht Bestandteil dieses Ansatzes. Nicht unterstützte Referenzen werden
sichtbar abgelehnt. Die Karten werden mit Build, Herkunft und Prüfsumme fest gebunden;
allein ein Dateiname reicht nicht als Identitätsnachweis.

Alternativen: Ein zusätzlicher hg19-Alignmentzweig würde Zeit, Speicher und eine weitere
Fehlerquelle hinzufügen. Eine reine Ergebnisimport-Funktion würde den gewünschten
automatischen Ablauf nicht erfüllen. Der direkte, buildgebundene Adapter ist deshalb
für das vorhandene GRCh38-Profil vorgesehen.

## Ablauf und Zustände

1. Bei Auswahl der Methylierung wird MARLIN ebenfalls angefordert. Der Desktop benennt
   dies ausdrücklich. Ohne Methylierung wird kein MARLIN-Modell geladen.
2. Vorprüfung: BAM-/Referenzidentität, unterstützte MM/ML-Information, installierte
   Modkit-Variante, freier Ausgabepfad, Modell-/Sonden-/Klassenidentität und verfügbare
   isolierte Laufzeit werden geprüft. Fehlende Voraussetzungen nennen konkrete Abhilfe.
3. Ein eigener MARLIN-Pileup verwendet die ausgewählten Modellsonden sowie kombinierte
   5mC/5hmC-Werte. Die vorhandenen regionalen 5mC-Werte werden nicht stillschweigend umgedeutet.
4. Pro Modellsonde wird `sum(N_mod) / sum(N_valid)` über die passenden Positionen gebildet.
   Die originalgetreue Merkmalsreihenfolge hat 357340 Einträge: beta >= 0.5 wird +1,
   beta < 0.5 wird -1, nicht beobachtet wird 0. Ein gemessener Nullwert bleibt beobachtet.
5. Die separate CPU-Laufzeit führt das Originalmodell ohne Training aus. Ein fester
   Threadplan, Ressourcen-/Laufzeitidentität und sämtliche Eingabe-/Ausgabeprüfsummen
   werden dokumentiert. Die Produktionsauswertung benötigt keine Netzwerkverbindung.
6. Kein beobachtetes Modellmerkmal führt zu NO_CALL, ohne Modellaufruf. Ein erfolgreiches
   Modell liefert genau 42 endliche Scores in der festgelegten Reihenfolge. Der bestehende
   Klassengruppierungsmechanismus und die feste 0.8-Schwelle werden nachvollziehbar angewandt.
   Unter der Schwelle bleibt die Zuordnung UNKNOWN; ein Score ist keine Diagnosewahrscheinlichkeit.
7. Fehlende Einrichtung wird NOT_RUN mit Ursache; Prozess-/Parserfehler werden FAILED.
   Andere auswertbare Module und der Bericht bleiben verfügbar. Ein MARLIN-Fehler ist
   kein negativer Befund und darf die übrige Genomauswertung nicht als fehlgeschlagen ausgeben.

## Datenvertrag, Wiederaufnahme und Bericht

Ein versionierter MARLIN-Laufbericht gehört zum bestehenden unveränderlichen Run-Verzeichnis.
Er enthält Run-/Probenidentität, Build, Status, Ursache, CpG-Zahl, Modell-/Laufzeitnachweise,
Scores, Entscheidung und den getrennten Validierungsstatus `UNVALIDATED_RESEARCH`.

Der Pipeline-Schritt und der Modulstatus werden explizit gespeichert. Nur ein aktueller,
prüfsummengeprüfter MARLIN-Lauf darf Ergebnisse liefern. Wiederaufnahme erfordert identische
Eingaben, Programmversion, Referenz, Sonden, Modell, Parameter und Laufzeitidentität.
Alte Dateien nach einem Fehler oder nach Abwahl bleiben als Diagnoseartefakte erhalten,
werden aber nicht in den neuen Bericht übernommen.

HTML mit und ohne JavaScript, JSON und XLSX zeigen dieselbe Entscheidung. Unterhalb der
Schwelle darf eine führende Klasse als Modellrangfolge erscheinen, jedoch nicht als
bestätigte Klassifikation. Fehlende Werte bleiben fehlend. Probenkennung, Modell und
Status müssen zwischen Ergebnisvertrag und MARLIN-Artefakt übereinstimmen.

Bereits archivierte Berichte werden nicht überschrieben. Ein nachträglicher MARLIN-Lauf
oder neu erzeugter Bericht erhält eine eigene Herkunftsangabe. Die bestehende
Forschungskennzeichnung und die manuelle Freigabe bleiben erhalten.

## Lokale Ressourcen und Bereitstellung

Vorhandene Ressourcen werden vor Übernahme erneut auf Byteidentität geprüft.
Quelle auf diesem Rechner: `Documents/ChatGPT/ONTSeq/work/extension-suite`.
Die folgenden Werte sind aus vorhandenen technischen Belegen übernommen und vor
Installation erneut zu verifizieren:

- Originalmodell, 1098245560 Bytes, SHA256
  `6210a674a0e7690b7b6c03184f5732a65037cf00ff4383e1892f26e90fdcb217`.
- Gepackte Laufzeit Python 3.10 / TensorFlow CPU 2.13.1 / Keras 2.13.1,
  SHA256 `b812927398645d3abaa3a56622bbf1b3279b70ff25067f42040d450083c822c6`.
- Merkmalsliste SHA256
  `9c271460d790d207ce91abaa28c6834d52dde8fe21dc4529033b14c7bb3cc4f5`.

Die Modelllaufzeit wird separat installiert; die übrigen Analysewerkzeuge erhalten keine
TensorFlow-Abhängigkeit. Aktive Analysen werden vor einer lokalen Aktualisierung erkannt
und nicht unterbrochen. Einstellungen, Ergebnisdateien und der bestehende Startpfad bleiben
verwendbar; die bisherige Installation wird für Rückkehr gesichert.

Offizielle Metadaten wurden lesend geprüft: Modell auf Zenodo 15565404, Version 1.0.0,
CC-BY-4.0, veröffentlichte MD5 `a12d4313ef7a97aa2df9776659bde7b2`; MARLIN-Code MIT.
Attributionen werden bei lokaler Paketierung mitgeführt. Das ist keine neue Lizenzierung
von ONTSeq. Es werden weder Patientendaten veröffentlicht noch Trainingsdaten benötigt.
Metadaten und historische Prüfberichte ersetzen keine erneute Prüfung der verwendeten Bytes.

## Technische Abnahme

- Parser: beide Builds, falsche Karte, doppelte Position/Strang, defekte MM/ML-Tags,
  ungleiche Tiefen, beta 0/0.49/0.5/1, fehlende Werte und leere Eingaben.
- Originalmodell: reale synthetische Ausführung in der vorgesehenen Laufzeit, Vergleich
  mit dem veröffentlichten Transformations-/Inferenzweg auf identischem Tensor;
  keine aus erfundenen Scores abgeleitete Freigabe.
- BAM bis Bericht: tatsächlich erzeugtes synthetisches BAM durch die installierte
  Modkit-Variante, gewichtete CpG-Werte, Tensor und Modellaufruf prüfen. Ein manuell
  geschriebenes bedMethyl allein ist kein Nachweis des vollständigen Übergangs.
- Zustände: NOT_RUN/FAILED/NO_CALL/UNKNOWN/erfolgreicher Modelllauf, falsche Probenidentität,
  veränderte Prüfsummen, Resume, Abwahl und veraltete Dateien.
- Oberfläche: Vorprüfung, automatischer Start, Fortschritt, Wiederanlauf und übereinstimmende
  statische/interaktive Ausgabe auf Windows/WSL; ohne Patientenanalyse.
- Repository: `make safety`, `make versions`, `make lint`, `make test`, relevante
  Windows-/Paketprüfungen und unabhängige Codeprüfung vor Integration in main.

Synthetische technische Abnahme belegt Ausführung und Datenübergabe. Eine analytische
Validierung an unabhängigen, fachlich bewerteten Proben bleibt davon getrennt und wird
weder aus der Installation noch aus den Softwaretests abgeleitet.

## Umsetzungseinheiten

1. Versionierte native Forschungsadapter-/Ressourcenverträge und buildgebundene Tests.
2. Kontrollierter Modkit-Pileup, Merkmalsbildung, isolierter Originalmodellaufruf und
   tatsächlicher synthetischer Übergangstest.
3. Pipeline, Signaturen, CLI/Service und Desktop-Vorprüfung mit automatischer Anforderung.
4. HTML/JSON/XLSX, echter Status und sichere Übernahme aktueller MARLIN-Ergebnisse.
5. Gesamtabnahme, Review, main-Integration und abgesicherte lokale Installation.

## Primärquellen und lokale technische Belege

- https://github.com/hovestadt/MARLIN/tree/442aa603415a54f62e7367794f9a31c6bc20fc2d
- https://github.com/hovestadt/MARLIN/blob/442aa603415a54f62e7367794f9a31c6bc20fc2d/MARLIN_realtime/2_process_pileup.R
- https://zenodo.org/records/15565404
- Lokale frühere Belege: `work/extension-suite/work/reports/marlin-probe.md` und
  `work/extension-suite/work/reports/marlin-bridge-acceptance.md` unter der genannten Quelle.
