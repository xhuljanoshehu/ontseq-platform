using System.Text.Json;
using System.Text.Json.Serialization;
using System.Windows.Media;

namespace ONTSeq.Desktop;

public sealed class DesktopSettings
{
    // Desktop installs must not require root privileges.  Core/CLI deliberately retains
    // /opt/ontseq as its independent server-style default.
    public const string DefaultResourceRootWsl =
        "~/.local/share/ontseq/resources-v" + DesktopVersion.CoreValue;

    public string WslDistribution { get; set; } = "Ubuntu";
    public string BackendCommand { get; set; } = "ontseq";
    public string? RuntimeBinWsl { get; set; }
    // Optional separately installed tool. The service checks its pinned version before
    // offering methylation; null preserves the backend's ordinary PATH lookup.
    public string? ModkitExecutableWsl { get; set; }
    public string ResourceRootWsl { get; set; } = DefaultResourceRootWsl;
    public string DefaultProfile { get; set; } = DesktopProfiles.DefaultProfileId;
    public string OutputDirectoryWindows { get; set; } = Path.Combine(
        Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments), "ONTSeq", "results");
    // Compatibility fields for one release. New profile runs resolve their resources from
    // ResourceRootWsl and never combine these explicit paths with a bundle context.
    public Dictionary<string, string> ReferenceLocksWsl { get; set; } = new(StringComparer.OrdinalIgnoreCase);
    public string? AdaptiveTargetBedWsl { get; set; }
    public string? AdaptiveTargetBedVersion { get; set; }
    public int Port { get; set; } = 8765;

    public bool HasAdaptiveTargetBedConfiguration =>
        !string.IsNullOrWhiteSpace(AdaptiveTargetBedWsl) ||
        !string.IsNullOrWhiteSpace(AdaptiveTargetBedVersion);

    public void ClearAdaptiveTargetBed()
    {
        AdaptiveTargetBedWsl = null;
        AdaptiveTargetBedVersion = null;
    }

    private static string DefaultUserSettingsPath => Path.Combine(
        Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
        "ONTSeq", "desktop.settings.json");

    public static string UserSettingsPath => SettingsOverridePath() ?? DefaultUserSettingsPath;

    private static string? SettingsOverridePath()
    {
        var value = Environment.GetEnvironmentVariable("ONTSEQ_DESKTOP_SETTINGS");
        if (string.IsNullOrWhiteSpace(value)) return null;
        var path = value.Trim();
        if (!Path.IsPathFullyQualified(path) || path.StartsWith(@"\\", StringComparison.Ordinal) ||
            path.StartsWith("//", StringComparison.Ordinal) || path.Any(char.IsControl) ||
            !string.Equals(Path.GetExtension(path), ".json", StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException(
                "ONTSEQ_DESKTOP_SETTINGS muss ein absoluter lokaler JSON-Dateipfad sein; " +
                "relative Pfade und UNC-Pfade sind nicht erlaubt.");
        return Path.GetFullPath(path);
    }

    public static DesktopSettings Load()
    {
        var settingsOverride = SettingsOverridePath();
        var candidates = settingsOverride is not null
            ? new[] { settingsOverride }
            : new[] {
            DefaultUserSettingsPath,
            Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.CommonApplicationData),
                "ONTSeq", "desktop.settings.json")
        };

        foreach (var file in candidates)
        {
            if (!File.Exists(file)) continue;
            var json = File.ReadAllText(file);
            var settings = JsonSerializer.Deserialize<DesktopSettings>(json, JsonDefaults.Options)
                           ?? throw new InvalidDataException($"Leere Desktop-Konfiguration: {file}");
            settings.ApplyProfileDefaults();
            return settings;
        }

        return new DesktopSettings();
    }

    public void ApplyProfileDefaults()
    {
        ResourceRootWsl = NormalizeResourceRootWsl(ResourceRootWsl);
        ModkitExecutableWsl = ValidateModkitExecutableWsl(ModkitExecutableWsl);
        ReferenceLocksWsl ??= new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        if (!DesktopProfiles.IsSupported(DefaultProfile))
            DefaultProfile = DesktopProfiles.DefaultProfileId;
    }

    public static string NormalizeResourceRootWsl(string? value)
    {
        var normalized = string.IsNullOrWhiteSpace(value)
            ? DefaultResourceRootWsl
            : value.Trim().Replace('\\', '/').TrimEnd('/');
        var isHomeRelative = normalized.StartsWith("~/", StringComparison.Ordinal);
        if ((!normalized.StartsWith("/", StringComparison.Ordinal) && !isHomeRelative) ||
            normalized is "/" or "~")
            throw new InvalidDataException(
                $"Der Resource-Root muss unter dem WSL-Home (z. B. {DefaultResourceRootWsl}) " +
                "oder als absoluter WSL-Pfad angegeben werden.");
        if (normalized.Any(char.IsControl) || normalized.Split('/').Any(part => part is "." or ".."))
            throw new InvalidDataException(
                "Der Resource-Root darf keine Steuerzeichen oder relativen Pfadsegmente enthalten.");
        return normalized;
    }

    public static string? ValidateModkitExecutableWsl(string? value)
    {
        if (value is null) return null;
        // Preserve the path, including spaces. It is one structured process argument,
        // never a command string, shell expression, or a Windows path to translate.
        if (!value.StartsWith("/", StringComparison.Ordinal) ||
            value.StartsWith("//", StringComparison.Ordinal) || value.EndsWith('/') ||
            value.Any(char.IsControl) || value.Split('/').Any(part => part is "." or "..") ||
            value.IndexOfAny(['\\', '$', '`', ';', '&', '|', '<', '>', '"', '\'', '(', ')']) >= 0)
            throw new InvalidDataException(
                "ModkitExecutableWsl muss null oder ein absoluter Linux-Dateipfad sein. " +
                "Verzeichniswurzeln, relative Segmente, Steuerzeichen und Shellausdrücke sind nicht erlaubt; " +
                "Leerzeichen innerhalb des Pfads bleiben erhalten.");
        return value;
    }

    public void SaveUserSettings()
    {
        ApplyProfileDefaults();
        var path = UserSettingsPath;
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);
        var temp = path + ".tmp";
        var json = JsonSerializer.Serialize(this, JsonDefaults.OptionsIndented);
        File.WriteAllText(temp, json + Environment.NewLine);
        File.Move(temp, path, overwrite: true);
    }

    public bool TryReferenceLockFor(string genomeBuild, out string referenceLock)
    {
        if (ReferenceLocksWsl.TryGetValue(genomeBuild, out var value) && !string.IsNullOrWhiteSpace(value))
        {
            referenceLock = value;
            return true;
        }
        referenceLock = string.Empty;
        return false;
    }

    public string ReferenceLockFor(string genomeBuild)
    {
        if (TryReferenceLockFor(genomeBuild, out var value)) return value;
        throw new InvalidOperationException(
            $"Für {genomeBuild} ist noch kein Reference-Lock konfiguriert. " +
            "Öffne 'System einrichten' und wähle die zum BAM passende FASTA/FAI-Referenz.");
    }
}

public sealed record DesktopAnalysisProfile(
    string ProfileId,
    string DisplayName,
    string GenomeBuild,
    string Assay,
    string DictionaryContract,
    bool AdaptiveSampling)
{
    public override string ToString() => DisplayName;

    public string DictionaryLabel =>
        string.Equals(DictionaryContract, "grch38_canonical_25", StringComparison.Ordinal)
            ? "Canonical-25 (chr1–22, chrX, chrY, chrM)"
            : string.Equals(
                DictionaryContract,
                "grch37_ucsc_hg19_canonical_25",
                StringComparison.Ordinal)
                ? "UCSC hg19 Canonical-25 (chr1–22, chrX, chrY, chrM=16571)"
            : GenomeBuild == "GRCh37"
                ? "vollständige GRCh37.p13 GENCODE 19 Assembly (inkl. Patches/Haplotypen)"
                : "vollständige GRCh38 Primary Assembly";
}

public static class DesktopProfiles
{
    public const string DefaultProfileId = "AML_LCWGS_GRCh38";
    public const string Grch37LcwgsProfileId = "AML_LCWGS_GRCh37";
    public const string Grch37UcscHg19Canonical25LcwgsProfileId =
        "AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25";
    public const string Grch37AdaptiveSamplingProfileId = "AML_AS_111_GRCh37";
    public const string Grch37UcscHg19Canonical25AdaptiveSamplingProfileId =
        "AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25";
    public const string AdaptiveSamplingProfileId = "AML_AS_111_GRCh38";
    public const string Canonical25LcwgsProfileId = "AML_LCWGS_GRCh38_CANONICAL25";
    public const string Canonical25AdaptiveSamplingProfileId =
        "AML_AS_111_GRCh38_CANONICAL25";

    public static readonly IReadOnlyList<DesktopAnalysisProfile> Supported =
    [
        new(
            DefaultProfileId,
            "AML low-coverage WGS · GRCh38",
            "GRCh38",
            "lcwgs",
            "exact_full",
            false),
        new(
            AdaptiveSamplingProfileId,
            "AML Adaptive Sampling 111 Gene · GRCh38",
            "GRCh38",
            "adaptive_sampling",
            "exact_full",
            true),
        new(
            Canonical25LcwgsProfileId,
            "AML low-coverage WGS · GRCh38 Canonical-25",
            "GRCh38",
            "lcwgs",
            "grch38_canonical_25",
            false),
        new(
            Canonical25AdaptiveSamplingProfileId,
            "AML Adaptive Sampling 111 Gene · GRCh38 Canonical-25",
            "GRCh38",
            "adaptive_sampling",
            "grch38_canonical_25",
            true),
        new(
            Grch37LcwgsProfileId,
            "AML low-coverage WGS · GRCh37.p13 / GENCODE 19",
            "GRCh37",
            "lcwgs",
            "exact_full",
            false),
        new(
            Grch37AdaptiveSamplingProfileId,
            "AML Adaptive Sampling · GRCh37.p13 (110/111 kartiert)",
            "GRCh37",
            "adaptive_sampling",
            "exact_full",
            true),
        new(
            Grch37UcscHg19Canonical25LcwgsProfileId,
            "AML low-coverage WGS · UCSC hg19 Canonical-25",
            "GRCh37",
            "lcwgs",
            "grch37_ucsc_hg19_canonical_25",
            false),
        new(
            Grch37UcscHg19Canonical25AdaptiveSamplingProfileId,
            "AML Adaptive Sampling · UCSC hg19 (110/111 kartiert)",
            "GRCh37",
            "adaptive_sampling",
            "grch37_ucsc_hg19_canonical_25",
            true)
    ];

    public static bool IsSupported(string? profileId) =>
        Supported.Any(profile => string.Equals(
            profile.ProfileId, profileId, StringComparison.Ordinal));

    public static DesktopAnalysisProfile Require(string profileId) =>
        Supported.SingleOrDefault(profile => string.Equals(
            profile.ProfileId, profileId, StringComparison.Ordinal))
        ?? throw new InvalidOperationException(
            $"Nicht unterstütztes Desktop-Profil: {profileId}. " +
            "Verfügbar sind buildgetrennte lcWGS- und Adaptive-Sampling-Profile für " +
            "GRCh38, natives GRCh37.p13 und UCSC hg19 Canonical-25.");
}

public static class JsonDefaults
{
    public static readonly JsonSerializerOptions Options = new()
    {
        PropertyNameCaseInsensitive = true,
        PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
        DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull
    };

    public static readonly JsonSerializerOptions OptionsIndented = new(Options)
    {
        WriteIndented = true
    };
}

public sealed record RunStartRequest(
    [property: JsonPropertyName("bam")] string Bam,
    [property: JsonPropertyName("sample_id")] string SampleId,
    [property: JsonPropertyName("run_id")] string? RunId,
    [property: JsonPropertyName("profile")] string Profile,
    [property: JsonPropertyName("genome_build")] string GenomeBuild,
    [property: JsonPropertyName("assay")] string Assay,
    [property: JsonPropertyName("target_bed")] string? TargetBed = null,
    [property: JsonPropertyName("target_bed_version")] string? TargetBedVersion = null,
    [property: JsonPropertyName("include_methylation")] bool IncludeMethylation = false);

public sealed record MethylationProbeRequest(
    [property: JsonPropertyName("bam_path")] string BamPath,
    [property: JsonPropertyName("force_refresh")] bool ForceRefresh = false);

public sealed record MethylationScanRequest(
    [property: JsonPropertyName("bam_path")] string BamPath);

public sealed record MethylationProbeResponse(
    [property: JsonPropertyName("bam_path")] string BamPath,
    [property: JsonPropertyName("status")] string Status,
    [property: JsonPropertyName("reason")] string Reason,
    [property: JsonPropertyName("checked_reads")] long CheckedReads,
    [property: JsonPropertyName("complete")] bool Complete,
    [property: JsonPropertyName("methylation_available")] bool MethylationAvailable = false,
    [property: JsonPropertyName("methylation_unavailable_reason")] string? MethylationUnavailableReason = null,
    [property: JsonPropertyName("expected_modkit_version")] string? ExpectedModkitVersion = null,
    [property: JsonPropertyName("reason_code")] string? ReasonCode = null,
    [property: JsonPropertyName("elapsed_seconds")] double ElapsedSeconds = 0,
    [property: JsonPropertyName("scan_mode")] string ScanMode = "quick",
    [property: JsonPropertyName("reader")] string? Reader = null,
    [property: JsonPropertyName("detected_modifications")] string[]? DetectedModifications = null)
{
    // Desktop-only context. A preparation error is not a zero-read scan result,
    // and this flag is never accepted from or sent to the service contract.
    [JsonIgnore]
    public bool PreparationFailed { get; init; }

    public MethylationProbeResponse RequireBam(string bamPath)
    {
        if (!string.Equals(BamPath, bamPath, StringComparison.Ordinal) ||
            Status is not ("detected" or "not_detected" or "unknown") ||
            string.IsNullOrWhiteSpace(Reason) || CheckedReads < 0 ||
            !double.IsFinite(ElapsedSeconds) || ElapsedSeconds < 0 ||
            ScanMode is not ("quick" or "thorough") ||
            (Status == "not_detected" && !Complete))
            throw new InvalidDataException(
                "Die Methylierungsprüfung hat keine eindeutig zur ausgewählten BAM passende Antwort geliefert.");
        return this;
    }
}

public static class MethylationProbeFailures
{
    public static MethylationProbeResponse FromException(
        string bamPath, Exception error, string scanMode, bool beforeProbe = false)
    {
        var preparationFailed = beforeProbe || error is WslPrerequisiteException;
        var code = error switch
        {
            WslPrerequisiteException prerequisite => prerequisite.ReasonCode,
            MethylationServiceException service => service.ReasonCode ?? "service_unavailable",
            System.Net.Http.HttpRequestException => "service_unavailable",
            TaskCanceledException or TimeoutException => "request_timeout",
            _ => preparationFailed ? "preparation_failed" : "worker_failed"
        };
        return new MethylationProbeResponse(bamPath, "unknown", error.Message, 0, false,
            ReasonCode: code, ScanMode: scanMode) { PreparationFailed = preparationFailed };
    }
}

public static class DesktopOutputDirectory
{
    public static void EnsureExists(string windowsPath)
    {
        try { _ = PathBridge.WindowsToWsl(windowsPath); }
        catch (Exception error) when (error is ArgumentException or InvalidOperationException or NotSupportedException)
        {
            throw new WslPrerequisiteException("output_path_unsupported",
                "Der konfigurierte Ausgabepfad muss ein vollständiger unterstützter Windows-Laufwerkspfad sein.\n" +
                $"Ausgabepfad (Windows): {windowsPath}\n{error.Message}",
                windowsPath: windowsPath, innerException: error);
        }
        try { Directory.CreateDirectory(windowsPath); }
        catch (Exception error) when (error is IOException or UnauthorizedAccessException or
            ArgumentException or NotSupportedException or System.Security.SecurityException)
        {
            throw new WslPrerequisiteException(
                error is ArgumentException or NotSupportedException ? "output_path_unsupported" : "output_unwritable",
                "Der Windows-Ausgabeordner konnte nicht vorbereitet werden. " +
                "Bitte einen nutzbaren Ordner prüfen oder auswählen.\n" +
                $"Ausgabeordner (Windows): {windowsPath}\n{error.Message}",
                windowsPath: windowsPath, innerException: error);
        }
    }
}

public sealed record MethylationScanSnapshot(
    [property: JsonPropertyName("scan_id")] string ScanId,
    [property: JsonPropertyName("bam_path")] string BamPath,
    [property: JsonPropertyName("state")] string State,
    [property: JsonPropertyName("checked_reads")] long CheckedReads,
    [property: JsonPropertyName("elapsed_seconds")] double ElapsedSeconds,
    [property: JsonPropertyName("result")] MethylationProbeResponse? Result)
{
    public static string RequireId(string scanId) =>
        scanId.Length == 32 && scanId.All(character => character is >= '0' and <= '9' or >= 'a' and <= 'f')
            ? scanId
            : throw new InvalidDataException("Ungültige Kennung der Methylierungsprüfung.");

    public MethylationScanSnapshot RequireScan(string bamPath, string? scanId = null)
    {
        RequireId(ScanId);
        if (!string.Equals(BamPath, bamPath, StringComparison.Ordinal) ||
            (scanId is not null && !string.Equals(ScanId, scanId, StringComparison.Ordinal)) ||
            State is not ("running" or "completed" or "cancelled" or "failed") ||
            CheckedReads < 0 || !double.IsFinite(ElapsedSeconds) || ElapsedSeconds < 0 ||
            (State == "running" && Result is not null) ||
            (State == "completed" && Result is null))
            throw new InvalidDataException("Der Prüfstatus passt nicht zur ausgewählten BAM oder Prüfung.");
        Result?.RequireBam(bamPath);
        return this;
    }
}

// A result is valid only for the exact selection event, profile and service instance.
// Selecting the same path again is a new event and never reuses an earlier decision.
public sealed record MethylationProbeScope(
    long Revision, string BamPath, string ProfileId, string? ServiceInstanceId);

public sealed class MethylationProbeGuard
{
    private long _revision;
    public void Invalidate() => _revision++;
    public MethylationProbeScope Capture(string bamPath, string profileId, string? serviceInstanceId) =>
        new(_revision, bamPath, profileId, serviceInstanceId);
    public bool IsCurrentSelection(MethylationProbeScope scope, string bamPath, string profileId) =>
        scope.Revision == _revision && scope.BamPath == bamPath && scope.ProfileId == profileId;
    public bool IsCurrent(MethylationProbeScope scope, string bamPath, string profileId, string? serviceInstanceId) =>
        IsCurrentSelection(scope, bamPath, profileId) && scope.ServiceInstanceId == serviceInstanceId;
}

public sealed record MethylationProbePresentation(string Title, string Explanation, string Progress, string Detail)
{
    public string ConfirmationText => $"{Title}\n\n{Explanation}\n\n{Progress}";

    public static string FormatProgress(long checkedReads, double elapsedSeconds) =>
        $"{checkedReads.ToString("N0", System.Globalization.CultureInfo.GetCultureInfo("de-DE"))} Reads geprüft · " +
        $"{elapsedSeconds.ToString("F1", System.Globalization.CultureInfo.GetCultureInfo("de-DE"))} Sekunden";

    public static MethylationProbePresentation From(MethylationProbeResponse probe)
    {
        if (probe.PreparationFailed) return PreparationFailure(probe);
        var title = probe.Status switch
        {
            "detected" => "5mC-Methylierungsinformationen erkannt",
            "not_detected" => "Keine Methylierungsinformationen erkannt",
            _ => probe.ReasonCode switch
            {
                "reader_unavailable" => "BAM-Lesewerkzeug nicht verfügbar",
                "input_unavailable" => "BAM nicht lesbar",
                "file_changed" => "BAM während der Prüfung verändert",
                "read_error" => "BAM-Lesefehler",
                "worker_failed" => "Prüfung fehlgeschlagen",
                "service_unavailable" => "Lokaler Dienst nicht erreichbar",
                "request_timeout" => "Dienst antwortet nicht rechtzeitig",
                "busy" => "Prüfung momentan nicht verfügbar",
                "cancelled" => "Prüfung abgebrochen",
                _ => "Methylierungsstatus offen"
            }
        };
        var explanation = probe.Status switch
        {
            "detected" => probe.MethylationAvailable
                ? "Die BAM enthält auswertbare 5mC-Tags. Vor jedem Analysestart entscheiden Sie ausdrücklich, ob Methylierung ergänzt wird."
                : "Die BAM enthält 5mC-Tags. Die zusätzliche Auswertung ist in dieser Installation noch nicht verfügbar. " +
                  (probe.MethylationUnavailableReason ?? "Werkzeug oder versionierte Policy sind nicht bestätigt."),
            "not_detected" => "Die gesamte BAM wurde geprüft; es wurden keine unterstützten Methylierungsinformationen gefunden.",
            _ => (probe.ReasonCode switch
            {
                "sample_incomplete" => "In der begrenzten Stichprobe wurden keine auswertbaren 5mC-Tags gefunden. Weitere Reads wurden noch nicht geprüft.",
                "timeout" => "Das Zeitlimit der Prüfung wurde erreicht, bevor ein abschließendes Ergebnis vorlag.",
                "cancelled" => "Die Prüfung wurde auf Ihren Wunsch beendet, bevor ein abschließendes Ergebnis vorlag.",
                "reader_unavailable" => "Das benötigte BAM-Lesewerkzeug steht im lokalen Dienst nicht zur Verfügung.",
                "input_unavailable" => "Die ausgewählte BAM kann vom lokalen Dienst nicht geöffnet werden.",
                "file_changed" => "Die Datei wurde während der Prüfung verändert. Das bisherige Prüfergebnis wurde verworfen.",
                "record_limit" => "Die Prüfgrenze wurde erreicht. Die BAM wurde noch nicht vollständig geprüft.",
                "invalid_tags" => "Es wurden unvollständige oder ungültige Modifikationstags gefunden.",
                "unsupported_modification" => "Es wurden Modifikationsinformationen gefunden, aber keine für diese Auswertung bestätigten 5mC-Tags.",
                "read_error" => "Beim Lesen der BAM trat ein Fehler auf. Die Prüfung konnte nicht abgeschlossen werden.",
                "worker_failed" => "Die Hintergrundprüfung wurde durch einen technischen Fehler beendet.",
                "service_unavailable" => "Die Verbindung zum lokalen Dienst ist fehlgeschlagen. Es liegt kein neues Prüfergebnis vor.",
                "request_timeout" => "Der lokale Dienst hat innerhalb des Zeitlimits keine Antwort auf die Prüfanfrage geliefert.",
                "busy" => "Der lokale Dienst ist bereits ausgelastet. Bitte nach Abschluss der laufenden Arbeit erneut prüfen.",
                _ => "Die Vorprüfung konnte kein abschließendes Ergebnis liefern."
            }) + " Damit ist keine sichere Aussage möglich, ob auswertbare Methylierungsinformationen vorhanden sind."
        };
        var detail = $"Prüfmodus: {(probe.ScanMode == "thorough" ? "gründlich" : "Stichprobe")} · " +
            $"Code: {probe.ReasonCode ?? "nicht angegeben"}\n{probe.Reason}";
        if (!string.IsNullOrWhiteSpace(probe.Reader)) detail += $"\nBAM-Leser: {probe.Reader}";
        return new(title, explanation, FormatProgress(probe.CheckedReads, probe.ElapsedSeconds), detail);
    }

    private static MethylationProbePresentation PreparationFailure(MethylationProbeResponse probe)
    {
        var (title, action) = probe.ReasonCode switch
        {
            "wsl_unavailable" => ("WSL nicht verfügbar", "In „System einrichten“ den WSL-Status prüfen."),
            "backend_unavailable" => ("ONTSeq-Laufzeit nicht verfügbar", "In „System einrichten“ die ONTSeq-Laufzeit prüfen."),
            "input_path_unsupported" => ("BAM-Pfad nicht unterstützt", "Die BAM über einen unterstützten Windows-Laufwerkspfad auswählen."),
            "output_path_unsupported" => ("Ausgabepfad nicht unterstützt", "Den konfigurierten Ausgabepfad auf einen vollständigen unterstützten Windows-Laufwerkspfad korrigieren."),
            "input_drive_unavailable" => ("BAM-Laufwerk in WSL nicht erreichbar", "Die Verbindung zum BAM-Laufwerk in WSL kann fehlen oder unterbrochen sein. Die Verfügbarkeit des Laufwerks und dessen Einbindung in WSL prüfen."),
            "input_root_missing" => ("BAM-Ordner in WSL nicht gefunden", "Den ausgewählten BAM-Ordner und seine Erreichbarkeit in WSL prüfen."),
            "input_root_unreadable" => ("BAM-Ordner in WSL nicht lesbar", "Die Leserechte und Zugriffsrechte für den BAM-Ordner in WSL prüfen."),
            "output_drive_unavailable" => ("Ausgabelaufwerk in WSL nicht erreichbar", "Die Verbindung zum Ausgabelaufwerk in WSL kann fehlen oder unterbrochen sein. Die Einbindung prüfen oder einen erreichbaren Ausgabeordner wählen."),
            "output_unwritable" => ("Ausgabeordner nicht nutzbar", "Einen beschreibbaren Ausgabeordner wählen. Der Pfad darf keine vorhandene Datei bezeichnen."),
            "resource_root_missing" => ("Ressourcenordner nicht gefunden", "In „System einrichten“ den konfigurierten Ressourcenordner und seinen Pfad prüfen."),
            "resource_root_unreadable" => ("Ressourcenordner nicht lesbar", "Die Zugriffsrechte für den konfigurierten Ressourcenordner prüfen."),
            "runtime_file_missing" => ("Benötigte Runtime-Datei fehlt", "In „System einrichten“ die ONTSeq-Laufzeit prüfen und bei Bedarf reparieren."),
            "runtime_file_unreadable" => ("Runtime-Datei nicht lesbar", "Die Zugriffsrechte der genannten Datei in der ONTSeq-Laufzeit prüfen."),
            "runtime_tool_unavailable" => ("Benötigtes Runtime-Werkzeug nicht verfügbar", "Das in den Prüfdetails genannte Werkzeug in der ONTSeq-Laufzeit prüfen."),
            "reference_lock_unavailable" => ("Referenzdefinition nicht verfügbar", "Die in den Prüfdetails genannte Referenzdefinition und ihre Zugriffsrechte prüfen."),
            "target_bed_unavailable" => ("Zielregionendatei nicht verfügbar", "Die für das Profil konfigurierte BED-Datei und ihre Zugriffsrechte prüfen."),
            "resource_profile_unavailable" => ("Analyseprofil nicht bereit", "In „System einrichten“ die Ressourcen des gewählten Referenzbuilds prüfen und die konkret gemeldete Lücke beheben."),
            "service_unavailable" => ("Lokaler Dienst nicht erreichbar", "Den lokalen Dienst und seine Verbindung prüfen, anschließend erneut prüfen."),
            "request_timeout" => ("Dienstvorbereitung nicht rechtzeitig bestätigt", "Den Status des lokalen Dienstes prüfen, bevor eine neue Prüfung angefordert wird."),
            "busy" => ("Lokaler Dienst bereits ausgelastet", "Die laufende Arbeit abschließen lassen, anschließend die BAM-Prüfung erneut anfordern."),
            "prerequisite_check_failed" => ("Vorbereitung nicht bestätigt", "Den konkreten Prüfschritt und die Pfade in den technischen Prüfdetails prüfen."),
            _ => ("BAM-Prüfung konnte nicht vorbereitet werden", "Den konkreten Fehler in den technischen Prüfdetails prüfen und die Konfiguration oder den lokalen Dienst korrigieren.")
        };
        return new(title,
            "Die BAM-Prüfung hat noch nicht begonnen. " + action +
            " Es ist noch keine Aussage über Methylierungsinformationen möglich.",
            "Prüfung nicht gestartet · noch keine BAM-Reads geprüft",
            $"Phase: Vorbereitung vor der BAM-Prüfung · Code: {probe.ReasonCode ?? "preparation_failed"}\n{probe.Reason}");
    }
}

public sealed record StageSnapshot(
    [property: JsonPropertyName("stage")] string Stage,
    [property: JsonPropertyName("title")] string Title,
    [property: JsonPropertyName("status")] string Status,
    [property: JsonPropertyName("reason")] string Reason,
    [property: JsonPropertyName("required")] bool Required,
    [property: JsonPropertyName("verification")] string Verification,
    [property: JsonPropertyName("resumed")] bool Resumed,
    [property: JsonPropertyName("duration_seconds")] double? DurationSeconds);

public sealed record RunJobResponse(
    [property: JsonPropertyName("run_id")] string RunId,
    [property: JsonPropertyName("sample_id")] string SampleId,
    [property: JsonPropertyName("state")] string State,
    [property: JsonPropertyName("detail")] string Detail,
    [property: JsonPropertyName("stages")] List<StageSnapshot> Stages,
    [property: JsonPropertyName("started_at")] string StartedAt,
    [property: JsonPropertyName("finished_at")] string? FinishedAt,
    [property: JsonPropertyName("profile")] string? Profile = null,
    [property: JsonPropertyName("detected_genome_build")] string? DetectedGenomeBuild = null);

public sealed record ServiceRootResponse(
    [property: JsonPropertyName("posix")] string Posix,
    [property: JsonPropertyName("display")] string Display);

public sealed record ServiceLaunchExpectation(
    string InstanceId,
    string ResourceRoot,
    string OutputDir,
    string AllowedRoot);

public sealed record ServiceConfigResponse(
    [property: JsonPropertyName("version")] string Version,
    [property: JsonPropertyName("output_dir")] string OutputDir,
    [property: JsonPropertyName("busy")] bool Busy,
    [property: JsonPropertyName("not_wired")] List<string> NotWired,
    [property: JsonPropertyName("profiles")] List<string>? Profiles = null,
    [property: JsonPropertyName("instance_id")] string? InstanceId = null,
    [property: JsonPropertyName("resource_root")] string? ResourceRoot = null,
    [property: JsonPropertyName("roots")] List<ServiceRootResponse>? Roots = null)
{
    public bool SupportsProfile(string profileId) =>
        Profiles?.Contains(profileId, StringComparer.Ordinal) == true;

    public ServiceConfigResponse RequireProfile(string profileId)
    {
        _ = DesktopProfiles.Require(profileId);
        if (SupportsProfile(profileId)) return this;
        var available = Profiles is { Count: > 0 }
            ? string.Join(", ", Profiles)
            : "keine";
        throw new InvalidOperationException(
            $"Der lokale Dienst kann das ausgewählte Profil {profileId} aus seinem aktuellen " +
            $"Resource-Root nicht auflösen. Vom Dienst gemeldete Profile: {available}. " +
            "Bitte 'System einrichten' öffnen und die betroffene Build-Familie reparieren.");
    }

    public ServiceConfigResponse RequireInstance(string expectedInstanceId)
    {
        if (!Guid.TryParseExact(expectedInstanceId, "N", out _))
            throw new ArgumentException(
                "Die erwartete Dienstinstanz muss als 32-stellige GUID ohne Trennzeichen angegeben werden.",
                nameof(expectedInstanceId));
        if (string.Equals(InstanceId, expectedInstanceId, StringComparison.Ordinal)) return this;
        throw new ServiceInstanceMismatchException(expectedInstanceId, InstanceId);
    }

    public ServiceConfigResponse RequireLaunch(ServiceLaunchExpectation expected)
    {
        RequireInstance(expected.InstanceId);
        var expectedResourceRoot = NormalizePosixPath(expected.ResourceRoot);
        var expectedOutputDir = NormalizePosixPath(expected.OutputDir);
        var expectedAllowedRoot = NormalizePosixPath(expected.AllowedRoot);
        var rootsMatch = Roots is { Count: 1 } &&
            string.Equals(
                NormalizePosixPath(Roots[0].Posix), expectedAllowedRoot, StringComparison.Ordinal);
        if (string.Equals(
                NormalizePosixPath(ResourceRoot), expectedResourceRoot, StringComparison.Ordinal) &&
            string.Equals(
                NormalizePosixPath(OutputDir), expectedOutputDir, StringComparison.Ordinal) &&
            rootsMatch)
            return this;
        throw new ServiceLaunchScopeMismatchException();
    }

    private static string NormalizePosixPath(string? value)
    {
        var normalized = value?.Trim() ?? "";
        return normalized.Length > 1 ? normalized.TrimEnd('/') : normalized;
    }
}

public sealed class ServiceInstanceMismatchException : InvalidOperationException
{
    public ServiceInstanceMismatchException(string expectedInstanceId, string? observedInstanceId)
        : base(
            "Am gewählten lokalen Port antwortet nicht die von diesem Desktop gestartete " +
            "ONTSeq-Instanz. Der fremde Dienst wird nicht übernommen oder beendet; ONTSeq " +
            "versucht einen neuen isolierten Port.")
    {
        ExpectedInstanceId = expectedInstanceId;
        ObservedInstanceId = observedInstanceId;
    }

    public string ExpectedInstanceId { get; }
    public string? ObservedInstanceId { get; }
}

public sealed class ServiceLaunchScopeMismatchException : InvalidOperationException
{
    public ServiceLaunchScopeMismatchException()
        : base(
            "Die vom lokalen Dienst gemeldeten Ressourcen-, Ausgabe- oder Eingabeverzeichnisse " +
            "stimmen nicht mit diesem Desktop-Start überein. Die Verbindung wird aus " +
            "Sicherheitsgründen nicht verwendet.")
    {
    }
}

public sealed record StageDisplay(string Title, string Status, string Reason)
{
    // Timeline colors mirror the portable report's status semantics: NO_CALL stays
    // neutral amber, NOT_RUN neutral gray; only FAILED is red, and a stage that has
    // not reported yet stays hollow.
    private static readonly Brush CompletedBrush = FrozenBrush("#2E7D5B");
    private static readonly Brush NoCallBrush = FrozenBrush("#B26A00");
    private static readonly Brush FailedBrush = FrozenBrush("#B3261E");
    private static readonly Brush NotRunBrush = FrozenBrush("#8A93A3");
    private static readonly Brush RunningBrush = FrozenBrush("#174A6E");
    private static readonly Brush PendingBrush = FrozenBrush("#9AA4B2");

    public string Symbol => Status switch
    {
        "COMPLETED" => "✓",
        "NO_CALL" => "○",
        "FAILED" => "✕",
        "NOT_RUN" => "—",
        "RUNNING" => "●",
        _ => "○"
    };

    public Brush DotFill => Status switch
    {
        "COMPLETED" => CompletedBrush,
        "NO_CALL" => NoCallBrush,
        "FAILED" => FailedBrush,
        "NOT_RUN" => NotRunBrush,
        "RUNNING" => RunningBrush,
        _ => Brushes.White
    };

    public Brush DotStroke => Status is "COMPLETED" or "NO_CALL" or "FAILED" or "NOT_RUN" or "RUNNING"
        ? DotFill
        : PendingBrush;

    public string StatusCaption => Status switch
    {
        "COMPLETED" => "abgeschlossen",
        "NO_CALL" => "kein verwertbarer Call",
        "FAILED" => "fehlgeschlagen",
        "NOT_RUN" => "nicht gelaufen",
        "RUNNING" => "läuft",
        _ => "ausstehend"
    };

    public string NodeToolTip => string.IsNullOrWhiteSpace(Reason)
        ? $"{Title} · {StatusCaption}"
        : $"{Title} · {StatusCaption}: {Reason}";

    private static Brush FrozenBrush(string hex)
    {
        var brush = (SolidColorBrush)new BrushConverter().ConvertFromInvariantString(hex)!;
        brush.Freeze();
        return brush;
    }
}
