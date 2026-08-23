using System.Text.Json;
using System.Text.Json.Serialization;

namespace ONTSeq.Desktop;

public sealed class DesktopSettings
{
    public string WslDistribution { get; set; } = "Ubuntu";
    public string BackendCommand { get; set; } = "ontseq";
    public string? RuntimeBinWsl { get; set; }
    public string OutputDirectoryWindows { get; set; } = Path.Combine(
        Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments), "ONTSeq", "results");
    public Dictionary<string, string> ReferenceLocksWsl { get; set; } = new(StringComparer.OrdinalIgnoreCase);
    public string? AdaptiveTargetBedWsl { get; set; }
    public string? AdaptiveTargetBedVersion { get; set; }
    public int Port { get; set; } = 8765;

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
            return JsonSerializer.Deserialize<DesktopSettings>(json, JsonDefaults.Options)
                   ?? throw new InvalidDataException($"Leere Desktop-Konfiguration: {file}");
        }

        return new DesktopSettings();
    }

    public void SaveUserSettings()
    {
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
    [property: JsonPropertyName("run_id")] string RunId,
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
    [property: JsonPropertyName("finished_at")] string? FinishedAt);

public sealed record ServiceConfigResponse(
    [property: JsonPropertyName("version")] string Version,
    [property: JsonPropertyName("output_dir")] string OutputDir,
    [property: JsonPropertyName("busy")] bool Busy,
    [property: JsonPropertyName("not_wired")] List<string> NotWired);

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
