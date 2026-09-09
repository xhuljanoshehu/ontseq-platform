using System.Net.Http;
using System.Net.Http.Json;
using System.Text.Json;
using System.Text.RegularExpressions;

namespace ONTSeq.Desktop;

public sealed class OntSeqServiceClient : IDisposable
{
    private const string TokenHeader = "X-ONTSeq-Token";
    private static readonly Regex TokenPattern = new(
        "const\\s+TOKEN\\s*=\\s*\"(?<token>[^\"]+)\"",
        RegexOptions.Compiled | RegexOptions.CultureInvariant);

    private readonly HttpClient _http;
    private bool _bootstrapped;

    public static Uri WorkspaceUri(int port, string? preferredProfileId = null)
    {
        var baseUri = $"http://127.0.0.1:{port}/workspace";
        if (string.IsNullOrWhiteSpace(preferredProfileId)) return new Uri(baseUri);
        var profile = DesktopProfiles.Require(preferredProfileId).ProfileId;
        return new Uri(baseUri + "?profile=" + Uri.EscapeDataString(profile));
    }

    public OntSeqServiceClient(int port)
    {
        _http = new HttpClient
        {
            BaseAddress = new Uri($"http://127.0.0.1:{port}/"),
            Timeout = TimeSpan.FromSeconds(20)
        };
    }

    public async Task<ServiceConfigResponse> BootstrapAsync(
        TimeSpan timeout,
        Func<bool>? backendExited,
        Func<string>? backendLog,
        ServiceLaunchExpectation expectation,
        CancellationToken cancellationToken)
    {
        _bootstrapped = false;
        _http.DefaultRequestHeaders.Remove(TokenHeader);
        var deadline = DateTimeOffset.UtcNow + timeout;
        Exception? last = null;
        while (DateTimeOffset.UtcNow < deadline)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (backendExited?.Invoke() == true)
                throw new InvalidOperationException("ONTSeq Backend wurde beim Start beendet.\n" + (backendLog?.Invoke() ?? ""));
            try
            {
                var html = await _http.GetStringAsync("", cancellationToken);
                var match = TokenPattern.Match(html);
                if (!match.Success) throw new InvalidDataException("Der ONTSeq-Dienst lieferte keinen Sitzungstoken.");
                _http.DefaultRequestHeaders.Remove(TokenHeader);
                _http.DefaultRequestHeaders.Add(TokenHeader, match.Groups["token"].Value);
                var config = (await GetJsonAsync<ServiceConfigResponse>(
                    "api/config", cancellationToken)).RequireLaunch(expectation);
                _bootstrapped = true;
                return config;
            }
            catch (Exception error) when (
                error is ServiceInstanceMismatchException or ServiceLaunchScopeMismatchException)
            {
                _bootstrapped = false;
                _http.DefaultRequestHeaders.Remove(TokenHeader);
                throw;
            }
            catch (Exception error) when (error is HttpRequestException or TaskCanceledException or InvalidDataException)
            {
                _bootstrapped = false;
                _http.DefaultRequestHeaders.Remove(TokenHeader);
                last = error;
                await Task.Delay(300, cancellationToken);
            }
        }

        throw new TimeoutException("ONTSeq Backend wurde nicht rechtzeitig erreichbar.", last);
    }

    public async Task<ServiceConfigResponse> GetConfigAsync(CancellationToken cancellationToken)
    {
        EnsureBootstrapped();
        return await GetJsonAsync<ServiceConfigResponse>("api/config", cancellationToken);
    }

    public async Task<RunJobResponse> StartRunAsync(RunStartRequest request, CancellationToken cancellationToken)
    {
        EnsureBootstrapped();
        using var response = await _http.PostAsJsonAsync("api/runs", request, JsonDefaults.Options, cancellationToken);
        var body = await response.Content.ReadAsStringAsync(cancellationToken);
        if (!response.IsSuccessStatusCode)
            throw new InvalidOperationException($"Analyse konnte nicht gestartet werden ({(int)response.StatusCode}): {ReadError(body)}");
        return JsonSerializer.Deserialize<RunJobResponse>(body, JsonDefaults.Options)
               ?? throw new InvalidDataException("Leere Antwort beim Start der Analyse.");
    }

    public async Task<MethylationProbeResponse> ProbeMethylationAsync(
        string bamPath, CancellationToken cancellationToken, bool forceRefresh = false)
    {
        EnsureBootstrapped();
        using var response = await _http.PostAsJsonAsync(
            "api/methylation/probe", new MethylationProbeRequest(bamPath, forceRefresh),
            JsonDefaults.Options, cancellationToken);
        var body = await response.Content.ReadAsStringAsync(cancellationToken);
        if (!response.IsSuccessStatusCode)
            throw MethylationError(response, body);
        return (JsonSerializer.Deserialize<MethylationProbeResponse>(body, JsonDefaults.Options)
            ?? throw new InvalidDataException("Leere Antwort bei der Methylierungsprüfung."))
            .RequireBam(bamPath);
    }

    public async Task<MethylationScanSnapshot> StartMethylationScanAsync(
        string bamPath, CancellationToken cancellationToken)
    {
        EnsureBootstrapped();
        using var response = await _http.PostAsJsonAsync(
            "api/methylation/scans", new MethylationScanRequest(bamPath),
            JsonDefaults.Options, cancellationToken);
        return (await ReadMethylationScanAsync(response, cancellationToken)).RequireScan(bamPath);
    }

    public async Task<MethylationScanSnapshot> GetMethylationScanAsync(
        string scanId, string bamPath, CancellationToken cancellationToken)
    {
        EnsureBootstrapped();
        using var response = await _http.GetAsync(
            "api/methylation/scans/" + MethylationScanSnapshot.RequireId(scanId), cancellationToken);
        return (await ReadMethylationScanAsync(response, cancellationToken)).RequireScan(bamPath, scanId);
    }

    public async Task<MethylationScanSnapshot> CancelMethylationScanAsync(
        string scanId, string bamPath, CancellationToken cancellationToken)
    {
        EnsureBootstrapped();
        using var response = await _http.PostAsJsonAsync(
            "api/methylation/scans/" + MethylationScanSnapshot.RequireId(scanId) + "/cancel",
            new { }, JsonDefaults.Options, cancellationToken);
        return (await ReadMethylationScanAsync(response, cancellationToken)).RequireScan(bamPath, scanId);
    }

    private static async Task<MethylationScanSnapshot> ReadMethylationScanAsync(
        HttpResponseMessage response, CancellationToken cancellationToken)
    {
        var body = await response.Content.ReadAsStringAsync(cancellationToken);
        if (!response.IsSuccessStatusCode) throw MethylationError(response, body);
        return JsonSerializer.Deserialize<MethylationScanSnapshot>(body, JsonDefaults.Options)
            ?? throw new InvalidDataException("Leere Antwort bei der gründlichen Methylierungsprüfung.");
    }

    private static MethylationServiceException MethylationError(HttpResponseMessage response, string body)
    {
        string? code = null;
        try
        {
            using var doc = JsonDocument.Parse(body);
            if (doc.RootElement.TryGetProperty("reason_code", out var value) && value.ValueKind == JsonValueKind.String)
                code = value.GetString();
        }
        catch (JsonException) { }
        return new MethylationServiceException(
            $"Methylierungsprüfung nicht verfügbar ({(int)response.StatusCode}): {ReadError(body)}", code);
    }

    public async Task<RunJobResponse> GetRunAsync(string runId, CancellationToken cancellationToken)
    {
        EnsureBootstrapped();
        return await GetJsonAsync<RunJobResponse>($"api/runs/{Uri.EscapeDataString(runId)}", cancellationToken);
    }

    private async Task<T> GetJsonAsync<T>(string relativeUri, CancellationToken cancellationToken)
    {
        using var response = await _http.GetAsync(relativeUri, cancellationToken);
        var body = await response.Content.ReadAsStringAsync(cancellationToken);
        if (!response.IsSuccessStatusCode)
            throw new InvalidOperationException($"ONTSeq Backend meldet {(int)response.StatusCode}: {ReadError(body)}");
        return JsonSerializer.Deserialize<T>(body, JsonDefaults.Options)
               ?? throw new InvalidDataException($"Leere Antwort von {relativeUri}.");
    }

    private static string ReadError(string body)
    {
        try
        {
            using var doc = JsonDocument.Parse(body);
            return doc.RootElement.TryGetProperty("error", out var error) ? error.GetString() ?? body : body;
        }
        catch (JsonException)
        {
            return body;
        }
    }

    private void EnsureBootstrapped()
    {
        if (!_bootstrapped) throw new InvalidOperationException("ONTSeq API ist noch nicht verbunden.");
    }

    public void Dispose() => _http.Dispose();
}

public sealed class MethylationServiceException(string message, string? reasonCode)
    : InvalidOperationException(message)
{
    public string? ReasonCode { get; } = reasonCode;
}

public static class ProvenanceReader
{
    public static async Task<string?> ReadGenomeBuildAsync(
        string outputRoot,
        string runId,
        string sampleId,
        CancellationToken cancellationToken)
    {
        var path = Path.Combine(outputRoot, runId, sampleId, "provenance", "run.json");
        if (!File.Exists(path)) return null;
        try
        {
            await using var stream = new FileStream(
                path,
                FileMode.Open,
                FileAccess.Read,
                FileShare.ReadWrite | FileShare.Delete);
            using var doc = await JsonDocument.ParseAsync(stream, cancellationToken: cancellationToken);
            return doc.RootElement.TryGetProperty("genome_build", out var build) &&
                   build.ValueKind == JsonValueKind.String
                ? build.GetString()
                : null;
        }
        catch (IOException)
        {
            return null;
        }
        catch (JsonException)
        {
            return null;
        }
    }

    public static async Task<List<StageSnapshot>> ReadStagesAsync(
        string outputRoot,
        string runId,
        string sampleId,
        CancellationToken cancellationToken)
    {
        var path = Path.Combine(outputRoot, runId, sampleId, "provenance", "run.json");
        if (!File.Exists(path)) return [];
        try
        {
            await using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete);
            using var doc = await JsonDocument.ParseAsync(stream, cancellationToken: cancellationToken);
            if (!doc.RootElement.TryGetProperty("stages", out var stages) || stages.ValueKind != JsonValueKind.Array)
                return [];
            var result = new List<StageSnapshot>();
            foreach (var stage in stages.EnumerateArray())
            {
                result.Add(new StageSnapshot(
                    stage.GetProperty("stage").GetString() ?? "unknown",
                    stage.GetProperty("title").GetString() ?? "Unbekannter Schritt",
                    stage.GetProperty("status").GetString() ?? "UNKNOWN",
                    stage.TryGetProperty("reason", out var reason) ? reason.GetString() ?? "" : "",
                    stage.TryGetProperty("required", out var required) && required.GetBoolean(),
                    stage.TryGetProperty("verification", out var verification) ? verification.GetString() ?? "" : "",
                    stage.TryGetProperty("resumed", out var resumed) && resumed.GetBoolean(),
                    stage.TryGetProperty("duration_seconds", out var duration) && duration.ValueKind == JsonValueKind.Number
                        ? duration.GetDouble()
                        : null));
            }
            return result;
        }
        catch (IOException)
        {
            return [];
        }
        catch (JsonException)
        {
            // Atomic replacement normally prevents partial JSON, but a network filesystem
            // may briefly expose a transition. The next poll will retry.
            return [];
        }
    }
}
