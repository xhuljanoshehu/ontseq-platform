# ONTSeq 0.7.0 — gemeinsame Plattform

Research Use Only. Lokaler Integrationskandidat; keine klinische Freigabe.

## Zusammengeführte Entwicklungsstände

- Gemeinsame Ausgangsbasis: `a66e709326ab0dac57af77407e6b3d6e58167e92`.
- Methylierungsentwicklung: `068216c37ec13f767ca79d8a8d3afe10b5c56379`.
- Plattformentwicklung: `d7189d9e732e825c011144d39647d699f5069d6a` plus der lokal
  vorhandene, noch nicht eingecheckte 0.6.2-Stand vom 8. September 2026.
- Integration auf `codex/ontseq-0.7-integration`. Die Originalordner wurden nicht
  überschrieben. Ein vorhandenes installiertes 0.6.2 wird nicht automatisch ersetzt.

Der 0.6.2-Arbeitsstand wurde vor der Zusammenführung mit Datei-Prüfsummen erfasst.
Die ursprünglichen Quellstände und das Zusammenführungsprotokoll liegen lokal im
ignorierten Verzeichnis `work/integration07/`. Die ausgelieferte Quellinventur wird
mit `scripts/build_local_metadata.py` erzeugt; ein Entwicklungsstand erhält keinen
erfundenen Release-Commit.

## Bedienablauf

1. BAM auswählen. Die Software prüft lokal und begrenzt auf MM/ML-Informationen.
2. Analyseprofil und Probenkennung prüfen.
3. Bei erkannten Methylierungsinformationen ausdrücklich entscheiden, ob diese
   zusätzlich ausgewertet werden sollen. Nach einem Dateitausch ist neu zu entscheiden.
4. Lauf starten und Fortschritt der angeforderten Module verfolgen.
5. Genomische Ergebnisse und gegebenenfalls regionale Methylierungswerte gemeinsam
   ansehen. JSON, HTML und XLSX behalten die Herkunft und fachlichen Einschränkungen.

Ein unvollständiger Scan bedeutet „noch nicht geklärt“. Fehlendes modkit ist eine
separate technische Voraussetzung und wird vor der gewünschten Zusatzauswertung
sichtbar. Die reine Genomanalyse kann weiterhin gestartet werden.

## Wissenschaftliche Grenzen

Die regionale Auswertung liefert deskriptive Methylierungsanteile. Mischungsanalyse,
Read-Holdout und unabhängige Validierung bleiben eigene Forschungsabläufe mit zwei
beziehungsweise vier Quellen. Ein einzelnes BAM liefert weder eine Kalibrationskohorte
noch eine validierte Tumor-, Blasten-, Zell- oder DNA-Massenfraktion.

Ergebnisschema 0.3.0, Methylierungsschemata, Adapterversionen und technische Policies
werden getrennt von der Produktversion 0.7.0 geführt. Die modkit-Pileup-Bindung bleibt
0.4.1, die davon unabhängige Einzelread-TSV-Schnittstelle bleibt 0.6.1.

## Laufzeit und Nachweise

Die Umgebungsvorgaben enthalten `ont-modkit=0.4.1`, passend zur vorhandenen Policy.
Der Paketname und die verfügbare Version sind in der offiziellen
[Bioconda-Rezeptdokumentation](https://bioconda.github.io/recipes/ont-modkit/README.html)
und im [ONT-Release](https://github.com/nanoporetech/modkit/releases/tag/v0.4.1)
nachprüfbar. Die jeweiligen Lizenzbedingungen bleiben anwendbar. Eine vorhandene
ältere Laufzeit wird dadurch nicht nachträglich verändert.

Auf dem Windows-Host fehlt `make`; die vier vorgeschriebenen Ziele werden versucht
und anschließend mit den unveränderten zugrunde liegenden Werkzeugen ausgeführt.
Reale Symlink-Tests bleiben auf unterstützten Hosts aktiv. Fehlt unter Windows die
erforderliche Berechtigung, wird genau dieser Umgebungsmangel als Skip ausgewiesen;
andere Dateisystemfehler werden nicht unterdrückt.

Der Windows-Prozessschutz verwendet jetzt eine nebenwirkungsfreie Handle-Prüfung.
Hintergrund: [Python dokumentiert die abweichende Windows-Signalbehandlung](https://docs.python.org/3/library/os.html#os.kill);
[Microsoft beschreibt die verwendete Prozess-Synchronisation](https://learn.microsoft.com/en-us/windows/win32/api/synchapi/nf-synchapi-waitforsingleobject).

Die abschließenden Software-, Paket- und Oberflächenprüfungen werden im lokalen
Integrationsprotokoll aufgezeichnet. Sie ersetzen keine analytische oder klinische
Validierung. Vorhandene Studienberichte werden nicht als 0.7-Evidenz umetikettiert.


## Abschlussprüfung vom 8. September 2026

- Windows: 1.233 pytest-Fälle gesammelt, Gesamtlauf erfolgreich (11 umgebungsabhängige
  Skips); das vorgeschriebene unittest-Ziel führt 1.154 Tests erfolgreich aus.
- Linux: 1.230 pytest-Tests und 347 Untertests bestanden; unittest mit 1.154 Tests
  erfolgreich. Der zunächst wegen PATH übersprungene reale R-Export wurde anschließend
  mit Rscript 4.4.3 geprüft (3 Tests, 4 Untertests bestanden). Verbleibend übersprungen:
  Windows-spezifische Prozessprüfung und ausdrücklich optionaler Vollreferenz-Smoke.
- Repository Safety, Versionsabgleich, Ruff, Formatierung und mypy bestanden.
- 26 Frontend-Vertragstests bestanden; das lokal gebaute HTML ist reproduzierbar und
  benötigt keine externen Assets. Browserprüfung bei 1440 × 980 und 390 × 844 Pixeln:
  Startseite, Ja/Nein, erneute Auswahl derselben BAM, fehlendes Werkzeug, vollständiger
  negativer Scan und unklare Tags geprüft. Keine horizontale Überlappung auf schmaler
  Ansicht. Der gefundene Fehler im leeren Anfangszustand wurde mit Regressionstest behoben.
- Eigenständiger Windows-x64-Build und Desktop-Vertragstests bestanden, ohne Warnungen.
- Echte synthetische CLI-Pipeline mit 27 künstlichen Reads und einer 250-kb-Testreferenz:
  samtools 1.24, cramino 1.3.0, Sniffles2 2.8.0, CuteSV 2.1.3 und modkit 0.4.1.
  QC meldet WARN; CNV und ISCN sind begründet NOT_RUN. Methylierung liefert 717 CpGs
  in zwei Regionen. 650 Sites erreichen die unveränderte Mindestabdeckung; die zweite
  Region bleibt mit zu geringer Abdeckung ohne Zahlenwert. JSON, HTML, XLSX und das
  unsignierte Prüfsummenpaket werden erzeugt. Ergebnis- und Methylierungs-API stimmen
  mit den typisierten Dateien überein; unautorisierter Zugriff wird abgewiesen.

Die Messwerte sind ausschließlich künstliche Software-Testfälle. Sie belegen weder eine
Tumorfraktion noch eine Vollgenom-, analytische oder klinische Validierung. Der Status
UNVERIFIED_ADAPTER wird durch einen lokalen Smoke-Test nicht aufgehoben. Die Policy-Notiz
wurde aktualisiert; Schwellenwerte und Rechenregeln wurden dabei nicht verändert.
Lokale Logs, Aufrufe und Hash-Prüfungen liegen in `work/integration07/`; diese ignorierten
Testausgaben gehören nicht zum Quellcode und enthalten keine Patientendaten.

## Lokaler Start und getrennte Abhängigkeiten

Eine nachfolgende Einrichtungskorrektur behebt falsch dekodierte Umlaute in
WSL-Diagnosen. Installations- und Reparaturaktionen zeigen den gewählten Build
ausdrücklich; eine anschließende Statusprüfung verdeckt die vorherige Fehlermeldung
nicht mehr. Der Inhalt lässt sich scrollen, sodass alle Einrichtungsaktionen auch in
kleineren Fenstern und bei langen Statusmeldungen erreichbar bleiben. GRCh37/hg19 und GRCh38/hg38
können im gleichen Resource-Root liegen,
bleiben aber durch ihre Profile, Bundles und Prüfsummen getrennt. Bereits vorhandene,
vollständig validierte lokale Bundles können einen fehlgeschlagenen Download ersetzen;
die DNS-Systemkonfiguration oder Referenz-Prüfsummen werden dafür nicht verändert.

`dist/ONTSeq-0.7.0-Windows-Testpaket/Start-ONTSeq-0.7.cmd` verwendet ausschließlich die
mitgelieferte `desktop.local.settings.json`. Damit erhält 0.7.0 einen eigenen Runtime-
und Ausgabeordner. Globale Desktop-Einstellungen von 0.6.2 werden nicht geladen.
Die vorhandenen hg19-Referenzressourcen werden lokal genutzt, ohne Build-Konvertierung.
Dieses lokal konfigurierte Paket setzt die auf diesem Rechner vorhandene WSL-Installation
und Referenzressourcen voraus; es ist kein universeller Installer mit Referenzdownload.

modkit 0.4.1 wurde separat vom offiziellen ONT-Release für den lokalen Forschungstest
bereitgestellt und über `modkitExecutableWsl` angebunden. Originalarchiv, Original-Lizenz,
Herkunft und Inhaltsinventur liegen getrennt unter `work/integration07/modkit041/`.
Das gemeinsame Basisruntime-Archiv wurde nicht um dieses Werkzeug erweitert. Die
[ONT Public License des Releases](https://github.com/nanoporetech/modkit/blob/v0.4.1/LICENCE.txt)
enthält eigene Nutzungs- und Weitergabebedingungen; dies ist keine Freigabe zur öffentlichen
Weitergabe des Gesamtpakets. Das Release bot keinen signierten Hersteller-Digest an;
lokal berechnete SHA256-Werte dokumentieren die erhaltenen Dateien, keine Signaturprüfung.
