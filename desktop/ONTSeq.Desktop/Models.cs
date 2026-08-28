using System.Text.Json;
using System.Text.Json.Serialization;

namespace ONTSeq.Desktop;

public sealed class DesktopSettings
{
    // Desktop installs must not require root privileges.  Core/CLI deliberately retains
    // /opt/ontseq as its independent server-style default.
    public const string DefaultResourceRootWsl = "~/.local/share/ontseq/resources";

    public string WslDistribution { get; set; } = "Ubuntu";
    public string BackendCommand { get; set; } = "ontseq";
    public string? RuntimeBinWsl { get; set; }
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

    public static string UserSettingsPath => Path.Combine(
        Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
        "ONTSeq", "desktop.settings.json");

    public static DesktopSettings Load()
    {
        var candidates = new[]
        {
            UserSettingsPath,
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
                "Der Resource-Root muss unter dem WSL-Home (z. B. ~/.local/share/ontseq/resources) " +
                "oder als absoluter WSL-Pfad angegeben werden.");
        if (normalized.Any(char.IsControl) || normalized.Split('/').Any(part => part is "." or ".."))
            throw new InvalidDataException(
                "Der Resource-Root darf keine Steuerzeichen oder relativen Pfadsegmente enthalten.");
        return normalized;
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
    bool AdaptiveSampling)
{
    public override string ToString() => DisplayName;
}

public static class DesktopProfiles
{
    public const string DefaultProfileId = "AML_LCWGS_GRCh38";
    public const string AdaptiveSamplingProfileId = "AML_AS_111_GRCh38";

    public static readonly IReadOnlyList<DesktopAnalysisProfile> Supported =
    [
        new(
            DefaultProfileId,
            "AML low-coverage WGS · GRCh38",
            "GRCh38",
            "lcwgs",
            false),
        new(
            AdaptiveSamplingProfileId,
            "AML Adaptive Sampling 111 Gene · GRCh38",
            "GRCh38",
            "adaptive_sampling",
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
            "Dieser Arbeitsstand aktiviert ausschließlich die getrennten GRCh38-Profile.");
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
    [property: JsonPropertyName("target_bed_version")] string? TargetBedVersion = null);

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

public sealed record ServiceConfigResponse(
    [property: JsonPropertyName("version")] string Version,
    [property: JsonPropertyName("output_dir")] string OutputDir,
    [property: JsonPropertyName("busy")] bool Busy,
    [property: JsonPropertyName("not_wired")] List<string> NotWired,
    [property: JsonPropertyName("profiles")] List<string>? Profiles = null);

public sealed record StageDisplay(string Title, string Status, string Reason)
{
    public string Symbol => Status switch
    {
        "COMPLETED" => "✓",
        "NO_CALL" => "○",
        "FAILED" => "✕",
        "NOT_RUN" => "—",
        "RUNNING" => "●",
        _ => "○"
    };
}
