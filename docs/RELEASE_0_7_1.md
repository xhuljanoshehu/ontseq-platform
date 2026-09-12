# ONTSeq 0.7.1 – Methylierungsprüfung

Lokaler Engineering-Kandidat vom 9. September 2026. Research Use Only; keine klinische Freigabe.

Neue Desktop-Einstellungen verwenden `~/.local/share/ontseq/resources-v0.7.1` als
Ressourcenordner. Der Name wird wie der Runtime-Präfix aus der Core-Version abgeleitet.
Gespeicherte Pfade bleiben erhalten; bestehende Installationen werden dadurch weder
verschoben noch erneut heruntergeladen. Ein bewusst umbenannter Ordner darf erst nach
abgeschlossenem Umzug und Prüfsummenprüfung als neuer Pfad ausgewählt werden; siehe
[Ressourcensystem](REFERENCE_SYSTEM.md#resource-root).

Die Versionsnummer im Ordner bezeichnet die Softwareinstallation. Wissenschaftliche
Identitäten bleiben getrennt: beispielsweise `GRCh37_GENCODE19_HG19_v2`,
`AML_AS_111_GRCh37_v1` und `HEMATOLOGY_GRCh37_v1`, einschließlich ihrer unveränderten
Prüfsummen. CLI-Standard `/opt/ontseq` und Ergebnisschema `0.3.0` bleiben erhalten.

## Bedienung

Nach der Auswahl einer BAM führt ONTSeq eine begrenzte Schnellprüfung aus. Wird unterstützte
Cytosin-Methylierung (5mC) erkannt, entscheiden Sie vor dem Start ausdrücklich, ob die
regionale Methylierung zusätzlich zur Genomanalyse ausgewertet werden soll.

Bei einer unvollständigen Stichprobe zeigt die Oberfläche den Grund und die geprüften
Reads an. **Gründlicher prüfen** startet eine sequenzielle Hintergrundprüfung. Sie zeigt
Fortschritt und Laufzeit und lässt sich abbrechen. Die Prüfung endet bei einem passenden
Fund vorzeitig; nur ein fehlerfreies Lesen bis zum Dateiende erlaubt einen vollständigen
Negativbefund. Andere Modifikationen, unvollständige Tags, Leseprobleme, Speichergrenzen
und ein Zeitlimit bleiben ausdrücklich unklar. Die gründliche Prüfung hat eine technische
Zeitgrenze von fünf Minuten und kann bei sehr großen Dateien unvollständig bleiben.

Ein BAM- oder Profilwechsel verwirft frühere Anzeigen und Entscheidungen und bricht eine
zugehörige Hintergrundprüfung ab. Während ein Scan oder seine Bereinigung läuft, wird
kein neuer Analyselauf begonnen. Die Vorstartprüfung liest die aktuelle Datei erneut;
ein zuvor gefundener Read wird anhand seiner internen Dateiposition frisch geprüft.

## Technik und Grenzen

Die korrigierte Desktop-Fassung meldet Fehler vor dem Prüfstart gesondert. WSL-Zugriff,
Eingabelaufwerk, Eingabeordner, Ausgabeordner, Ressourcen und Runtime-Komponenten werden
mit eigenen Fehlercodes und passenden Handlungshinweisen geprüft. Eine fehlende oder
unterbrochene Laufwerkseinbindung ist kein Fehler des BAM-Lesers und kein Anlass, ein
Referenzbundle erneut zu installieren. Die Anwendung bindet Laufwerke nicht automatisch
neu ein. Nach einer gezielten Wiederherstellung kann die Prüfung erneut gestartet werden.

Der isolierte Reader verwendet pysam 0.24.0 direkt auf BAM-Datensätzen. Die Schnellprüfung
kann BAI-/CSI-Indizes für mehrere Stichprobenpositionen verwenden; eine beschädigte oder
nicht verfügbare Indexstichprobe darf keine Abwesenheit behaupten. Die bisherigen
SAM-Textzeilen werden nicht mehr erzeugt. Der Worker ist zeitlich und im Speicher begrenzt.
Die unterstützten Analyseumgebungen deklarieren pysam explizit; ein reines Core-Paket
ohne diesen optionalen Reader meldet dessen Fehlen verständlich.

Probe-Provenienz: `mm-ml-presence-v2`, konkreter Grundcode, Scanmodus, Reader, Laufzeit,
Read-Zahl, Vollständigkeit und Datei-Metadatenfingerprint. Dieser Fingerprint ist kein
Datei-Inhaltschecksum. Read-Namen, Sequenzen und interne Trefferpositionen werden nicht
über die Prüf-API ausgegeben. Anzeige-Caches sind begrenzt; sie ersetzen keine neue
Bestätigung vor dem Start oder die regulären Eingabeprüfungen der Analyse.

GRCh37/hg19 und GRCh38/hg38 bleiben getrennte Referenzfamilien. Die Änderung betrifft
Erkennung und Bedienung. Methylierungsrechenregeln, technische Schwellen, Referenzbundles
und das Ergebnisschema 0.3.0 bleiben unverändert. Die erforderliche
Validierungsfolgenbewertung steht in `CLINICAL_VALIDATION.md`.

## Prüfung

Synthetische Tests umfassen späte Treffer, lange Reads, fehlende/unvollständige Tags,
andere Modifikationen, vollständigen EOF, BAI/CSI, Dateiwechsel, Worker-Abbruch und
Speichergrenzen. Service- und Desktop-Prüfungen decken Authentifizierung, Pfadgrenzen,
Cache-Verfall, alte Antworten und bestätigten Abbruch ab. Die Browserprüfung umfasst
Desktop- und schmale Ansichten sowie Auswahlwechsel während einer verzögerten Antwort.
Reale Patienten-BAMs werden dafür weder gesucht noch geöffnet.
