using System.Diagnostics;
using System.Net;
using System.Net.Sockets;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;

namespace ONTSeq.Desktop;

public sealed record LoopbackPortSelection(int Port, bool UsedFallback);

public sealed class WslServiceLauncher : IAsyncDisposable
{
    private const string BaseLinuxPath = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin";
    private const string ReleaseVersion = DesktopVersion.CoreValue;

    private static class RuntimeAssets
    {
        public const string QcPolicy = "configs/qc/defaults.yaml";
        public const string TargetCoveragePolicy = "configs/qc/adaptive_target_coverage.technical.yaml";
        public const string Components = "configs/components/default.yaml";
        public const string SnifflesPolicy = "configs/sv/sniffles2.conservative.technical.yaml";
        public const string CuteSvPolicy = "configs/sv/cutesv.conservative.technical.yaml";
        public const string SvConsensusPolicy = "configs/sv/sniffles2_cutesv.consensus.technical.yaml";
        public const string SvEvidencePolicy = "configs/sv/evidence-priority.technical.yaml";
        public const string CnvPolicy = "configs/cnv/qdnaseq_ace.technical.yaml";
        public const string QdnaSeqScript = "scripts/run_qdnaseq_ace.R";

        public static IReadOnlyList<string> RequiredFiles { get; } = Array.AsReadOnly(
            new[]
            {
                QcPolicy,
                TargetCoveragePolicy,
                Components,
                SnifflesPolicy,
                CuteSvPolicy,
                SvConsensusPolicy,
                SvEvidencePolicy,
                CnvPolicy,
                QdnaSeqScript
            });

        public static IReadOnlyList<string> RequiredTools { get; } = Array.AsReadOnly(
            new[] { "ontseq", "Rscript", "samtools", "cramino", "sniffles", "cuteSV", "mosdepth" });
    }

    public const string Grch38ReferenceBundleId = "GRCh38_GENCODE50_MANE1.5_v1";
    public const string Grch37ReferenceBundleId = "GRCh37_GENCODE19_HG19_v2";
    public const string Grch37HematologyKnowledgeBundleId = "HEMATOLOGY_GRCh37_v1";
    public const string HematologyKnowledgeBundleId = "HEMATOLOGY_v3";
    public const string AmlAdaptivePanelBundleId = "AML_AS_111_GRCh38_v1";
    public const string Grch37AmlAdaptivePanelBundleId = "AML_AS_111_GRCh37_v1";
    public static IReadOnlyList<string> ManagedGrch38ResourceBundleIds { get; } =
        Array.AsReadOnly(new[] {
            Grch38ReferenceBundleId,
            HematologyKnowledgeBundleId,
            AmlAdaptivePanelBundleId
        });
    public static IReadOnlyList<string> ManagedGrch37ResourceBundleIds { get; } =
        Array.AsReadOnly(new[] {
            Grch37ReferenceBundleId,
            Grch37HematologyKnowledgeBundleId,
            Grch37AmlAdaptivePanelBundleId
        });

    public static IReadOnlyList<string> ManagedResourceBundleIds(string genomeBuild) =>
        genomeBuild switch
        {
            "GRCh38" => ManagedGrch38ResourceBundleIds,
            "GRCh37" => ManagedGrch37ResourceBundleIds,
            _ => throw new ArgumentOutOfRangeException(nameof(genomeBuild), genomeBuild,
                "Nicht unterstützter Genome-Build.")
        };
    private Process? _process;
    private Task? _stdoutPump;
    private Task? _stderrPump;
    private readonly object _outputLock = new();
    private readonly StringBuilder _stderr = new();
    private readonly StringBuilder _stdout = new();

    public string DiagnosticLog
    {
        get { lock (_outputLock) return $"{_stdout}\n{_stderr}".Trim(); }
    }

    public static LoopbackPortSelection SelectLoopbackPort(int preferredPort)
    {
        if (preferredPort is < 1 or > 65535)
            throw new ArgumentOutOfRangeException(
                nameof(preferredPort), preferredPort, "Der lokale Dienstport muss zwischen 1 und 65535 liegen.");

        if (CanBindLoopback(preferredPort))
            return new LoopbackPortSelection(preferredPort, false);

        // Ask Windows for an unused ephemeral IPv4 loopback port. The ONTSeq service binds
        // 127.0.0.1 as well, so an existing Desktop/WSL service is observed here. Releasing
        // the probe immediately before process start avoids taking over or stopping it.
        return SelectEphemeralLoopbackPort();
    }

    internal static LoopbackPortSelection SelectEphemeralLoopbackPort()
    {
        using var listener = new TcpListener(IPAddress.Loopback, 0);
        listener.Server.ExclusiveAddressUse = true;
        listener.Start();
        var selected = ((IPEndPoint)listener.LocalEndpoint).Port;
        return new LoopbackPortSelection(selected, true);
    }

    internal static string CreateServiceInstanceId() => Guid.NewGuid().ToString("N");

    internal static bool IsServiceLaunchCollision(Exception error, string? diagnosticLog = null)
    {
        if (error is ServiceInstanceMismatchException) return true;
        var detail = error.Message + "\n" + diagnosticLog;
        return new[]
        {
            "address already in use",
            "errno 98",
            "10048",
            "only one usage of each socket address",
            "failed to bind"
        }.Any(marker => detail.Contains(marker, StringComparison.OrdinalIgnoreCase));
    }

    private static bool CanBindLoopback(int port)
    {
        try
        {
            using var listener = new TcpListener(IPAddress.Loopback, port);
            listener.Server.ExclusiveAddressUse = true;
            listener.Start();
            return true;
        }
        catch (SocketException)
        {
            return false;
        }
    }

    public static string WorkspaceAllowedRootWindows(string? bam, string outputDirectory)
    {
        var root = !string.IsNullOrWhiteSpace(bam) && File.Exists(bam) &&
                   bam.EndsWith(".bam", StringComparison.OrdinalIgnoreCase)
            ? Path.GetDirectoryName(Path.GetFullPath(bam))!
            : Path.GetFullPath(outputDirectory);
        _ = PathBridge.WindowsToWsl(root);
        return root;
    }

    public static bool IsWithinAllowedRoot(string candidate, string? allowedRoot)
    {
        if (allowedRoot is null) return false;
        var root = Path.GetFullPath(allowedRoot).TrimEnd(Path.DirectorySeparatorChar);
        var path = Path.GetFullPath(candidate).TrimEnd(Path.DirectorySeparatorChar);
        return path.Equals(root, StringComparison.OrdinalIgnoreCase) ||
            path.StartsWith(root + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase);
    }

    public async Task VerifyPrerequisitesAsync(
        DesktopSettings settings,
        string allowedRootWindows,
        string referenceLockWsl,
        string expectedGenomeBuild,
        CancellationToken cancellationToken)
    {
        await VerifyWslBackendAsync(settings, cancellationToken);
        var checks = WorkspacePrerequisiteChecks(settings, allowedRootWindows).ToList();
        checks.Add(new WslPrerequisiteCheck(
            "reference_lock_unavailable", referenceLockWsl, null,
            "Die konfigurierte Referenzdatei ist in WSL nicht lesbar. " +
            "Prüfe die ausgewählte Referenzdatei in 'System einrichten'.\n" +
            $"Referenzdatei: {referenceLockWsl}",
            $"test -f {ShellQuote(referenceLockWsl)} && test -r {ShellQuote(referenceLockWsl)}"));
        AddRuntimePrerequisiteChecks(settings, checks);
        if (!string.IsNullOrWhiteSpace(settings.AdaptiveTargetBedWsl))
            checks.Add(new WslPrerequisiteCheck(
                "target_bed_unavailable", settings.AdaptiveTargetBedWsl, null,
                "Die konfigurierte Zielregionsdatei fehlt, ist leer oder nicht lesbar. " +
                $"Prüfe die ausgewählte BED-Datei.\nZielregionen: {settings.AdaptiveTargetBedWsl}",
                $"test -s {ShellQuote(settings.AdaptiveTargetBedWsl)} && test -r {ShellQuote(settings.AdaptiveTargetBedWsl)}"));
        await VerifyPathPrerequisitesAsync(settings, checks, cancellationToken);

        var reference = await CheckReferenceAsync(
            settings, referenceLockWsl, expectedGenomeBuild, cancellationToken);
        if (!reference.Ok)
            throw new WslPrerequisiteException("reference_lock_unavailable",
                $"Die gespeicherte {expectedGenomeBuild}-Referenz ist ungültig oder " +
                "unvollständig. Öffne 'System einrichten' und wähle den vollständigen " +
                $"FAI-Index der BAM-Referenz.\n{reference.Detail}", referenceLockWsl);
    }

    public async Task<(bool Ok, string Detail)> CheckWslAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken)
    {
        try
        {
            var status = await RunWslAsync(
                settings.WslDistribution, ["sh", "-lc", "printf ready"], cancellationToken);
            if (status.ExitCode == 0 && status.StdOut.Contains("ready", StringComparison.Ordinal))
                return (true, $"WSL-Distribution '{settings.WslDistribution}' ist erreichbar.");
            return (false, $"WSL2 bzw. die Distribution '{settings.WslDistribution}' ist nicht einsatzbereit. {status.StdErr}".Trim());
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested) { throw; }
        catch (Exception error)
        {
            return (false, "WSL2 konnte nicht gestartet werden: " + error.Message);
        }
    }

    public async Task<(bool Ok, string Detail)> CheckBackendAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken)
    {
        try
        {
            var versionProbe = await RunWslAsync(
                settings.WslDistribution,
                BackendInvocation(settings, "--version"),
                cancellationToken);
            if (versionProbe.ExitCode != 0 ||
                !string.Equals(versionProbe.StdOut.Trim(), ReleaseVersion, StringComparison.Ordinal))
            {
                var observed = string.IsNullOrWhiteSpace(versionProbe.StdOut)
                    ? versionProbe.StdErr.Trim()
                    : versionProbe.StdOut.Trim();
                return (
                    false,
                    $"Installierte ONTSeq Runtime '{observed}' entspricht nicht v{ReleaseVersion}. " +
                    "Bitte 'Runtime installieren' erneut ausführen.");
            }

            var result = await RunWslAsync(
                settings.WslDistribution,
                BackendInvocation(settings, "--help"),
                cancellationToken);
            if (result.ExitCode != 0)
                return (false, "ONTSeq Backend ist in WSL nicht einsatzbereit. " + result.StdErr);

            var referenceCapability = await RunWslAsync(
                settings.WslDistribution,
                BackendInvocation(settings, "validate-reference", "--help"),
                cancellationToken);
            if (referenceCapability.ExitCode != 0)
                return (
                    false,
                    "Die installierte ONTSeq Runtime ist veraltet. Bitte 'Runtime " +
                    $"installieren' ausführen, um sie auf Desktop/Core v{ReleaseVersion} zu aktualisieren.");

            var serviceCapability = await RunWslAsync(
                settings.WslDistribution,
                BackendInvocation(settings, "serve", "--help"),
                cancellationToken);
            if (serviceCapability.ExitCode != 0 ||
                !serviceCapability.StdOut.Contains("--target-coverage-policy", StringComparison.Ordinal) ||
                !serviceCapability.StdOut.Contains("--components", StringComparison.Ordinal) ||
                !serviceCapability.StdOut.Contains("--cutesv-policy", StringComparison.Ordinal) ||
                !serviceCapability.StdOut.Contains("--sv-consensus-policy", StringComparison.Ordinal) ||
                !serviceCapability.StdOut.Contains("--sv-evidence-policy", StringComparison.Ordinal) ||
                !serviceCapability.StdOut.Contains("--resource-root", StringComparison.Ordinal) ||
                !serviceCapability.StdOut.Contains("--instance-id", StringComparison.Ordinal))
            {
                return (
                    false,
                    $"Die installierte ONTSeq Runtime enthält nicht den vollständigen v{ReleaseVersion}-" +
                    "Desktop-Vertrag für Target Coverage, Komponentenauswahl und SV-Policies. Bitte " +
                    "'Runtime installieren' erneut ausführen.");
            }

            var resourceCapability = await RunWslAsync(
                settings.WslDistribution,
                BackendInvocation(settings, "references", "--help"),
                cancellationToken);
            if (resourceCapability.ExitCode != 0)
            {
                return (
                    false,
                    "Die installierte ONTSeq Runtime unterstützt noch keine manifestierten " +
                    "Referenz-Bundles. Bitte 'Runtime installieren' erneut ausführen.");
            }

            return (true, $"ONTSeq Backend v{ReleaseVersion} gefunden: {settings.BackendCommand}");
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested) { throw; }
        catch (Exception error)
        {
            return (false, "ONTSeq Backend fehlt oder kann nicht gestartet werden: " + error.Message);
        }
    }

    public async Task<(bool Ok, string Detail)> CheckReferenceAsync(
        DesktopSettings settings,
        string referenceLockWsl,
        string expectedGenomeBuild,
        CancellationToken cancellationToken)
    {
        var result = await RunWslAsync(
            settings.WslDistribution,
            BackendInvocation(
                settings,
                "validate-reference",
                referenceLockWsl,
                "--expected-genome-build",
                expectedGenomeBuild,
                "--require-canonical-assembly"),
            cancellationToken);
        if (result.ExitCode == 0 && !string.IsNullOrWhiteSpace(result.StdOut))
            return (true, result.StdOut.Trim());

        var detail = string.IsNullOrWhiteSpace(result.StdErr) ? result.StdOut : result.StdErr;
        return (false, detail.Trim());
    }

    public async Task<IReadOnlyDictionary<string, ResourceFamilyState>> CheckResourceFamiliesAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken)
    {
        settings.ApplyProfileDefaults();
        var args = ResourceManagementArguments(
            "status", settings.ResourceRootWsl).ToList();
        args.Add("--json");
        var result = await RunWslAsync(
            settings.WslDistribution,
            BackendInvocation(settings, args.ToArray()),
            cancellationToken);
        return InterpretResourceFamilyStates(result.ExitCode, result.StdOut, result.StdErr);
    }

    public async Task<(bool Ok, string Detail)> CheckResourceBundlesAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken,
        string genomeBuild = "GRCh38")
    {
        _ = ManagedResourceBundleIds(genomeBuild);
        var families = await CheckResourceFamiliesAsync(settings, cancellationToken);
        var family = families[genomeBuild];
        return (family.CanAnalyze, family.Detail);
    }

    public async Task<(bool Ok, string Detail)> ValidateResourceBundlesAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken) =>
        await RunResourceQueryAsync(settings, "validate", cancellationToken);

    private async Task<(bool Ok, string Detail)> RunResourceQueryAsync(
        DesktopSettings settings,
        string action,
        CancellationToken cancellationToken)
    {
        settings.ApplyProfileDefaults();
        var result = await RunWslAsync(
            settings.WslDistribution,
            BackendInvocation(
                settings,
                ResourceManagementArguments(
                    action, settings.ResourceRootWsl).ToArray()),
            cancellationToken);
        var detail = string.IsNullOrWhiteSpace(result.StdOut) ? result.StdErr : result.StdOut;
        return (result.ExitCode == 0, detail.Trim());
    }

    public static (bool Ok, string Detail) InterpretResourceStatus(
        int exitCode,
        string stdout,
        string stderr,
        string genomeBuild = "GRCh38")
    {
        _ = ManagedResourceBundleIds(genomeBuild);
        var family = InterpretResourceFamilyStates(exitCode, stdout, stderr)[genomeBuild];
        return (family.CanAnalyze, family.Detail);
    }

    public static IReadOnlyDictionary<string, ResourceFamilyState> InterpretResourceFamilyStates(
        int exitCode,
        string stdout,
        string stderr)
    {
        if (exitCode != 0)
        {
            var failure = (string.IsNullOrWhiteSpace(stderr) ? stdout : stderr).Trim();
            return DesktopResourcePolicy.UnavailableFamilies(
                string.IsNullOrWhiteSpace(failure)
                    ? "'ontseq references status' ist fehlgeschlagen."
                    : "Ressourcenstatus nicht verfügbar: " + failure);
        }

        try
        {
            using var document = JsonDocument.Parse(stdout);
            var root = document.RootElement;
            if (root.ValueKind != JsonValueKind.Object ||
                !root.TryGetProperty("references", out var references) ||
                references.ValueKind != JsonValueKind.Array ||
                !root.TryGetProperty("profiles", out var profileArray) ||
                profileArray.ValueKind != JsonValueKind.Array)
            {
                throw new InvalidDataException(
                    "Statusantwort benötigt die Arrays 'references' und 'profiles'.");
            }

            var installedReferences = new Dictionary<string, (bool Valid, string? GenomeBuild)>(
                StringComparer.Ordinal);
            foreach (var item in references.EnumerateArray())
            {
                if (item.ValueKind != JsonValueKind.Object ||
                    !item.TryGetProperty("bundle_id", out var id) ||
                    id.ValueKind != JsonValueKind.String ||
                    id.GetString() is not { } bundleId ||
                    !item.TryGetProperty("valid", out var valid) ||
                    valid.ValueKind is not (JsonValueKind.True or JsonValueKind.False))
                {
                    throw new InvalidDataException("Ungültiger Referenzeintrag in der Statusantwort.");
                }
                string? declaredBuild = null;
                if (item.TryGetProperty("genome_build", out var build))
                {
                    if (build.ValueKind != JsonValueKind.String)
                        throw new InvalidDataException("Ungültiger Genome-Build in der Statusantwort.");
                    declaredBuild = build.GetString();
                }
                if (!installedReferences.TryAdd(
                        bundleId,
                        (valid.ValueKind == JsonValueKind.True, declaredBuild)))
                    throw new InvalidDataException($"Doppelter Referenzeintrag: {bundleId}.");
            }

            var readyProfiles = new HashSet<string>(StringComparer.Ordinal);
            foreach (var item in profileArray.EnumerateArray())
            {
                if (item.ValueKind != JsonValueKind.String || item.GetString() is not { } profileId)
                    throw new InvalidDataException("Ungültiger Profileintrag in der Statusantwort.");
                readyProfiles.Add(profileId);
            }

            var observedProfiles = new HashSet<string>(readyProfiles, StringComparer.Ordinal);
            if (root.TryGetProperty("profile_status", out var profileStatus))
            {
                if (profileStatus.ValueKind != JsonValueKind.Array)
                    throw new InvalidDataException("'profile_status' muss ein Array sein.");
                foreach (var item in profileStatus.EnumerateArray())
                {
                    if (item.ValueKind != JsonValueKind.Object ||
                        !item.TryGetProperty("profile_id", out var id) ||
                        id.ValueKind != JsonValueKind.String || id.GetString() is not { } profileId)
                    {
                        throw new InvalidDataException(
                            "Ungültiger Profile-Status-Eintrag in der Statusantwort.");
                    }
                    observedProfiles.Add(profileId);
                }
            }

            var families = new Dictionary<string, ResourceFamilyState>(StringComparer.Ordinal);
            foreach (var genomeBuild in DesktopResourcePolicy.GenomeBuilds)
            {
                var managedBundles = ManagedResourceBundleIds(genomeBuild);
                var referenceBundleId = managedBundles[0];
                var requiredProfiles = DesktopProfiles.Supported
                    .Where(profile => profile.GenomeBuild == genomeBuild)
                    .Select(profile => profile.ProfileId)
                    .ToArray();
                var referencePresent = installedReferences.TryGetValue(
                    referenceBundleId, out var reference);
                var familyProfilesObserved = requiredProfiles.Any(observedProfiles.Contains);
                if (!referencePresent && !familyProfilesObserved)
                {
                    families[genomeBuild] = new ResourceFamilyState(
                        genomeBuild,
                        ResourceFamilyAvailability.NotInstalled,
                        $"Nicht installiert: {referenceBundleId} und die zugehörigen Profile fehlen. " +
                        "Die Installation benötigt Internetzugang und mehrere GB Speicherplatz.");
                    continue;
                }

                var issues = new List<string>();
                if (!referencePresent)
                    issues.Add(referenceBundleId + " fehlt");
                else if (!reference.Valid)
                    issues.Add(referenceBundleId + " ist ungültig");
                else if (reference.GenomeBuild is not null && !string.Equals(
                             reference.GenomeBuild, genomeBuild, StringComparison.Ordinal))
                    issues.Add(referenceBundleId + " deklariert " + reference.GenomeBuild);

                issues.AddRange(requiredProfiles
                    .Where(profile => !readyProfiles.Contains(profile))
                    .Select(profile => profile + " fehlt oder ist ungültig"));
                if (issues.Count > 0)
                {
                    families[genomeBuild] = new ResourceFamilyState(
                        genomeBuild,
                        ResourceFamilyAvailability.Incomplete,
                        "Unvollständig: " + string.Join(", ", issues) +
                        ". Reparatur erforderlich; fehlende Referenzdateien können einen Download benötigen.");
                    continue;
                }

                families[genomeBuild] = new ResourceFamilyState(
                    genomeBuild,
                    ResourceFamilyAvailability.Ready,
                    "Lokal vorhanden; Analyse ohne Download möglich (Schnellprüfung). Ressourcen: " +
                    string.Join(", ", managedBundles) + ". Profile: " +
                    string.Join(", ", requiredProfiles) + ".");
            }
            return families;
        }
        catch (Exception error) when (error is JsonException or InvalidDataException)
        {
            return DesktopResourcePolicy.UnavailableFamilies(
                "Ungültige Antwort von 'ontseq references status': " + error.Message);
        }
    }

    public async Task<string> InstallGrch38ProfileResourcesAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken) =>
        await RunResourceManagementAsync(
            settings,
            "install",
            Grch38ReferenceBundleId,
            cancellationToken);

    public async Task<string> RepairGrch38ProfileResourcesAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken) =>
        await RunResourceManagementAsync(
            settings,
            "repair",
            Grch38ReferenceBundleId,
            cancellationToken);

    public Task<string> InstallProfileResourcesAsync(
        DesktopSettings settings, string genomeBuild, CancellationToken cancellationToken) =>
        RunResourceManagementAsync(settings, "install", ManagedResourceBundleIds(genomeBuild)[0],
            cancellationToken);

    public Task<string> RepairProfileResourcesAsync(
        DesktopSettings settings, string genomeBuild, CancellationToken cancellationToken) =>
        RunResourceManagementAsync(settings, "repair", ManagedResourceBundleIds(genomeBuild)[0],
            cancellationToken);

    private async Task<string> RunResourceManagementAsync(
        DesktopSettings settings,
        string action,
        string bundleId,
        CancellationToken cancellationToken)
    {
        settings.ApplyProfileDefaults();
        var result = await RunWslAsync(
            settings.WslDistribution,
            BackendInvocation(
                settings,
                ResourceManagementArguments(
                    action, settings.ResourceRootWsl, bundleId).ToArray()),
            cancellationToken);
        if (result.ExitCode != 0)
        {
            var detail = string.IsNullOrWhiteSpace(result.StdErr) ? result.StdOut : result.StdErr;
            throw new InvalidOperationException(
                $"Resource-Bundle konnte nicht mit '{action}' verarbeitet werden.\n{detail.Trim()}");
        }
        return (string.IsNullOrWhiteSpace(result.StdOut) ? result.StdErr : result.StdOut).Trim();
    }

    public static IReadOnlyList<string> ResourceManagementArguments(
        string action,
        string resourceRootWsl,
        string? bundleId = null)
    {
        if (action is not ("status" or "list" or "validate" or "install" or "repair" or "import"))
            throw new ArgumentOutOfRangeException(nameof(action), action, "Unbekannter Resource-Befehl.");
        var requiresTarget = action is "install" or "repair" or "import";
        if (requiresTarget != !string.IsNullOrWhiteSpace(bundleId))
            throw new ArgumentException(
                requiresTarget
                    ? $"'{action}' benötigt eine Bundle-ID bzw. einen Importpfad."
                    : $"'{action}' akzeptiert kein Bundle-Ziel.",
                nameof(bundleId));

        var args = new List<string> { "references", action };
        if (requiresTarget) args.Add(bundleId!);
        args.Add("--resource-root");
        args.Add(DesktopSettings.NormalizeResourceRootWsl(resourceRootWsl));
        return args;
    }

    public async Task VerifyProfilePrerequisitesAsync(
        DesktopSettings settings,
        string allowedRootWindows,
        string profileId,
        CancellationToken cancellationToken)
    {
        var profile = DesktopProfiles.Require(profileId);
        settings.ApplyProfileDefaults();

        await VerifyWslBackendAsync(settings, cancellationToken);
        await VerifyPathPrerequisitesAsync(
            settings, ProfilePrerequisiteChecks(settings, allowedRootWindows), cancellationToken);

        IReadOnlyDictionary<string, ResourceFamilyState> families;
        try { families = await CheckResourceFamiliesAsync(settings, cancellationToken); }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested) { throw; }
        catch (Exception error)
        {
            throw new WslPrerequisiteException("wsl_unavailable",
                $"Der Ressourcenstatus konnte in WSL '{settings.WslDistribution}' nicht abgefragt werden. " +
                "Prüfe die Erreichbarkeit der Distribution und starte die Vorprüfung erneut.", innerException: error);
        }
        var family = families[profile.GenomeBuild];
        if (family.Availability == ResourceFamilyAvailability.Unavailable)
            throw new WslPrerequisiteException("prerequisite_check_failed",
                $"Der Status der {profile.GenomeBuild}-Ressourcen konnte nicht sicher bestimmt werden. " +
                "Prüfe zunächst die WSL-/Backend-Verbindung und frage den Ressourcenstatus erneut ab.\n" +
                family.Detail, settings.ResourceRootWsl);
        if (!family.CanAnalyze)
            throw new WslPrerequisiteException("resource_profile_unavailable",
                $"Die manifestierten {profile.GenomeBuild}-Profile sind nicht vollständig einsatzbereit. " +
                "Nutze in 'System einrichten' die vollständige Ressourcen-Reparatur für " +
                string.Join(", ", ManagedResourceBundleIds(profile.GenomeBuild)) + ".\n" +
                family.Detail, settings.ResourceRootWsl);
    }

    private async Task VerifyWslBackendAsync(DesktopSettings settings, CancellationToken cancellationToken)
    {
        var wsl = await CheckWslAsync(settings, cancellationToken);
        if (!wsl.Ok)
            throw new WslPrerequisiteException("wsl_unavailable", wsl.Detail);
        var backend = await CheckBackendAsync(settings, cancellationToken);
        if (!backend.Ok)
            throw new WslPrerequisiteException("backend_unavailable",
                backend.Detail + $"\nBackend: {settings.BackendCommand}\nRuntime-Werkzeuge: {settings.RuntimeBinWsl ?? "PATH"}",
                settings.BackendCommand);
    }

    private static string MappedWorkspacePath(string windowsPath, bool output)
    {
        try { return PathBridge.WindowsToWsl(windowsPath); }
        catch (Exception error) when (error is ArgumentException or InvalidOperationException or NotSupportedException)
        {
            throw new WslPrerequisiteException(
                output ? "output_path_unsupported" : "input_path_unsupported",
                (output ? "Der Ausgabepfad" : "Der ausgewählte BAM-Ordner") +
                " kann nicht nach WSL übersetzt werden. Wähle einen vollständigen Windows-Laufwerkspfad " +
                "(z. B. E:\\Daten). UNC- und bereits übersetzte Linux-Pfade werden hier nicht automatisch zugeordnet.",
                windowsPath: windowsPath, innerException: error);
        }
    }

    private static IReadOnlyList<WslPrerequisiteCheck> WorkspacePrerequisiteChecks(
        DesktopSettings settings, string allowedRootWindows)
    {
        var root = MappedWorkspacePath(allowedRootWindows, output: false);
        var output = MappedWorkspacePath(settings.OutputDirectoryWindows, output: true);
        var inputDrive = root[..6]; // PathBridge guarantees /mnt/<drive>, without inspecting files.
        var outputDrive = output[..6];
        var distribution = settings.WslDistribution;
        return new[]
        {
            new WslPrerequisiteCheck("input_drive_unavailable", inputDrive, allowedRootWindows,
                $"Das Laufwerk des ausgewählten BAM-Ordners ist in WSL '{distribution}' nicht erreichbar. " +
                "Prüfe, ob das Laufwerk in Windows verbunden und in dieser WSL-Distribution zugänglich ist. " +
                "Auch eine getrennte oder veraltete Einbindung kann diesen Fehler verursachen. " +
                "Bei einer abweichenden WSL-Laufwerkszuordnung muss der Pfad zuerst passend eingerichtet werden.\n" +
                $"Windows-Ordner: {allowedRootWindows}\nErwartetes WSL-Laufwerk: {inputDrive}",
                $"test -d {ShellQuote(inputDrive)}"),
            new WslPrerequisiteCheck("input_root_missing", root, allowedRootWindows,
                $"Der ausgewählte BAM-Ordner ist in WSL '{distribution}' nicht erreichbar. " +
                "Prüfe den Ordnerpfad und die Verbindung des Laufwerks; wähle den Ordner danach erneut.\n" +
                $"Windows-Ordner: {allowedRootWindows}\nErwarteter WSL-Ordner: {root}",
                $"test -d {ShellQuote(root)}"),
            new WslPrerequisiteCheck("input_root_unreadable", root, allowedRootWindows,
                "Der ausgewählte BAM-Ordner kann in WSL nicht gelesen oder betreten werden. " +
                $"Prüfe die Zugriffsrechte des Laufwerks und Ordners.\nBAM-Ordner in WSL: {root}",
                $"test -r {ShellQuote(root)} && test -x {ShellQuote(root)}"),
            new WslPrerequisiteCheck("output_drive_unavailable", outputDrive, settings.OutputDirectoryWindows,
                $"Das Ausgabelaufwerk ist in WSL '{distribution}' nicht erreichbar. " +
                "Prüfe seine Verbindung und WSL-Zuordnung oder wähle einen erreichbaren Ausgabeordner.\n" +
                $"Windows-Ausgabe: {settings.OutputDirectoryWindows}\nErwartetes WSL-Laufwerk: {outputDrive}",
                $"test -d {ShellQuote(outputDrive)}"),
            new WslPrerequisiteCheck("output_unwritable", output, settings.OutputDirectoryWindows,
                "Der Ausgabeordner kann in WSL nicht angelegt oder beschrieben werden. " +
                "Prüfe freien Speicher, Zugriffsrechte und ob am Pfad bereits eine Datei liegt; " +
                $"wähle bei Bedarf einen anderen Ausgabeordner.\nAusgabeordner in WSL: {output}",
                $"mkdir -p {ShellQuote(output)} && test -d {ShellQuote(output)} && test -w {ShellQuote(output)} && test -x {ShellQuote(output)}")
        };
    }

    internal static IReadOnlyList<WslPrerequisiteCheck> ProfilePrerequisiteChecks(
        DesktopSettings settings, string allowedRootWindows)
    {
        var checks = WorkspacePrerequisiteChecks(settings, allowedRootWindows).ToList();
        var resource = settings.ResourceRootWsl;
        checks.Add(new WslPrerequisiteCheck("resource_root_missing", resource, null,
            "Der konfigurierte Ressourcenordner ist in WSL nicht erreichbar. " +
            "Prüfe in 'System einrichten' zunächst den Resource-Root und dessen Laufwerk. " +
            "Vorhandene Ressourcen können dort ausgewählt werden.\n" + $"Resource-Root: {resource}",
            $"test -d {ShellPathExpression(resource)}"));
        checks.Add(new WslPrerequisiteCheck("resource_root_unreadable", resource, null,
            "Der Ressourcenordner kann in WSL nicht gelesen oder betreten werden. " +
            $"Prüfe seine Zugriffsrechte.\nResource-Root: {resource}",
            $"test -r {ShellPathExpression(resource)} && test -x {ShellPathExpression(resource)}"));
        AddRuntimePrerequisiteChecks(settings, checks);
        return checks;
    }

    private static void AddRuntimePrerequisiteChecks(DesktopSettings settings, List<WslPrerequisiteCheck> checks)
    {
        if (string.IsNullOrWhiteSpace(settings.RuntimeBinWsl)) return;
        foreach (var path in RequiredRuntimeFiles(settings))
        {
            checks.Add(new WslPrerequisiteCheck("runtime_file_missing", path, null,
                "Eine benötigte Runtime-Datei fehlt. Wähle in 'System einrichten' die aktuelle " +
                $"ONTSeq-Runtime oder installiere das lokale Runtime-Paket erneut.\nRuntime-Datei: {path}",
                $"test -f {ShellQuote(path)}"));
            checks.Add(new WslPrerequisiteCheck("runtime_file_unreadable", path, null,
                "Eine benötigte Runtime-Datei ist nicht lesbar. " +
                $"Prüfe die Zugriffsrechte der ausgewählten Runtime.\nRuntime-Datei: {path}",
                $"test -r {ShellQuote(path)}"));
        }
        foreach (var path in RequiredRuntimeTools(settings))
            checks.Add(new WslPrerequisiteCheck("runtime_tool_unavailable", path, null,
                "Ein benötigtes Runtime-Werkzeug fehlt oder ist nicht ausführbar. " +
                "Prüfe die ausgewählte Runtime und ihre Zugriffsrechte; bei unvollständiger Runtime " +
                $"kann das lokale Runtime-Paket erneut installiert werden.\nWerkzeug: {path}",
                $"test -f {ShellQuote(path)} && test -x {ShellQuote(path)}"));
    }

    internal static string PrerequisiteCommand(IReadOnlyList<WslPrerequisiteCheck> checks)
    {
        // Only a fixed numeric check identifier crosses stdout. No file contents or shell
        // diagnostics are read or echoed; one failing dependency stops subsequent actions.
        var commands = checks.Select((check, index) =>
            $"if ! ( {check.TestCommand} ) 2>/dev/null; then printf 'ONTSEQ_PREREQUISITE:%s\\n' '{index}'; exit 42; fi");
        return string.Join("\n", commands) + "\nprintf '%s\\n' 'ONTSEQ_PREREQUISITES_OK'";
    }

    internal static void RequirePrerequisites(
        IReadOnlyList<WslPrerequisiteCheck> checks, int exitCode, string stdout)
    {
        var output = stdout.Trim();
        if (exitCode == 0 && output == "ONTSEQ_PREREQUISITES_OK") return;
        const string prefix = "ONTSEQ_PREREQUISITE:";
        if (exitCode == 42 && output.StartsWith(prefix, StringComparison.Ordinal) &&
            int.TryParse(output[prefix.Length..], out var index) && index >= 0 && index < checks.Count)
            throw checks[index].Failure();
        throw new WslPrerequisiteException("prerequisite_check_failed",
            "Die WSL-Vorprüfung lieferte keinen eindeutigen Status. Prüfe, ob die ausgewählte " +
            "WSL-Distribution weiterhin erreichbar ist, und starte die Prüfung erneut.");
    }

    private static async Task VerifyPathPrerequisitesAsync(
        DesktopSettings settings, IReadOnlyList<WslPrerequisiteCheck> checks, CancellationToken cancellationToken)
    {
        try
        {
            var check = await RunWslAsync(settings.WslDistribution,
                ["sh", "-lc", PrerequisiteCommand(checks)], cancellationToken);
            RequirePrerequisites(checks, check.ExitCode, check.StdOut);
        }
        catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested) { throw; }
        catch (WslPrerequisiteException) { throw; }
        catch (Exception error)
        {
            throw new WslPrerequisiteException("wsl_unavailable",
                $"Die Pfadprüfung konnte in WSL '{settings.WslDistribution}' nicht ausgeführt werden. " +
                "Prüfe die Erreichbarkeit dieser Distribution und starte die Prüfung erneut.", innerException: error);
        }
    }

    public async Task<string> InstallBundledRuntimeAsync(
        DesktopSettings settings,
        string runtimeArchiveWindows,
        CancellationToken cancellationToken)
    {
        if (!File.Exists(runtimeArchiveWindows))
            throw new FileNotFoundException(
                "Das gebündelte Linux-Runtime-Paket fehlt. Bitte den vollständigen ONTSeq-Desktop-ZIP entpacken, nicht nur die EXE.",
                runtimeArchiveWindows);

        // The base conda-pack archive remains unchanged. Verify the matching Core wheel
        // separately before creating a new prefix; never modify the installed old runtime.
        var package = await RuntimePackage.VerifyAsync(runtimeArchiveWindows, ReleaseVersion, cancellationToken);

        var wsl = await CheckWslAsync(settings, cancellationToken);
        if (!wsl.Ok) throw new InvalidOperationException(wsl.Detail);

        var homeResult = await RunWslAsync(
            settings.WslDistribution, ["sh", "-lc", "printf %s \"$HOME\""], cancellationToken);
        if (homeResult.ExitCode != 0 || string.IsNullOrWhiteSpace(homeResult.StdOut))
            throw new InvalidOperationException("WSL-Home-Verzeichnis konnte nicht bestimmt werden. " + homeResult.StdErr);

        var home = homeResult.StdOut.Trim();
        var target = home + $"/.local/share/ontseq/runtime-v{ReleaseVersion}-{Guid.NewGuid():N}";
        var bin = target + "/bin";
        var runtimePath = bin + ":" + BaseLinuxPath;
        var installedRuntime = new DesktopSettings { RuntimeBinWsl = bin };
        var command =
            RuntimeInstallCoreCommand(package, target) + " && " +
            BundledRuntimePrerequisiteCommand(installedRuntime) + " && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/Rscript")} -e " +
            ShellQuote("stopifnot(requireNamespace('QDNAseq',quietly=TRUE), requireNamespace('QDNAseq.hg19',quietly=TRUE), requireNamespace('QDNAseq.hg38',quietly=TRUE), requireNamespace('ACE',quietly=TRUE))") +
            " && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} validate-reference --help >/dev/null && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} references --help >/dev/null && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--target-coverage-policy' && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--components' && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--cutesv-policy' && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--sv-consensus-policy' && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--sv-evidence-policy' && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--resource-root' && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--instance-id'";
        var install = await RunWslAsync(
            settings.WslDistribution, ["sh", "-lc", command], cancellationToken);
        if (install.ExitCode != 0)
            throw new InvalidOperationException(
                $"ONTSeq Linux-Runtime v{ReleaseVersion} konnte nicht installiert werden.\n" + install.StdErr);

        settings.RuntimeBinWsl = bin;
        settings.BackendCommand = bin + "/ontseq";
        settings.SaveUserSettings();
        return target;
    }

    public static string RuntimeInstallCoreCommand(RuntimePackage package, string target)
    {
        var normalized = DesktopSettings.NormalizeResourceRootWsl(target);
        if (!normalized.StartsWith('/') ||
            !PosixDirectoryName(normalized).EndsWith("/.local/share/ontseq", StringComparison.Ordinal) ||
            !normalized[(normalized.LastIndexOf('/') + 1)..].StartsWith(
                $"runtime-v{ReleaseVersion}-", StringComparison.Ordinal))
            throw new InvalidDataException("Die Runtime benötigt einen neuen, versionsgebundenen Installationspfad.");
        var bin = normalized + "/bin";
        var runtimePath = bin + ":" + BaseLinuxPath;
        return $"test ! -e {ShellQuote(normalized)} && mkdir -p {ShellQuote(normalized)} && " +
            $"tar -xzf {ShellQuote(PathBridge.WindowsToWsl(package.ArchivePath))} -C {ShellQuote(normalized)} && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/conda-unpack")} && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/python")} -I -m pip --isolated install " +
            "--no-deps --no-index --force-reinstall --no-cache-dir --disable-pip-version-check " +
            ShellQuote(PathBridge.WindowsToWsl(package.WheelPath)) + " && " +
            RuntimeVersionCheckCommand(normalized);
    }

    public static string RuntimeVersionCheckCommand(string target)
    {
        var bin = DesktopSettings.NormalizeResourceRootWsl(target) + "/bin";
        var script =
            "import subprocess,sys; " +
            "value=subprocess.check_output([sys.prefix+'/bin/ontseq','--version'],text=True).strip(); " +
            $"assert value == '{ReleaseVersion}', value; print(value)";
        return $"env PATH={ShellQuote(bin + ":" + BaseLinuxPath)} " +
            $"{ShellQuote(bin + "/python")} -I -c {ShellQuote(script)}";
    }

    public async Task<string> ConfigureReferenceAsync(
        DesktopSettings settings,
        string sourceWindows,
        string genomeBuild,
        CancellationToken cancellationToken)
    {
        if (!File.Exists(sourceWindows))
            throw new FileNotFoundException("Referenzdatei nicht gefunden.", sourceWindows);

        var backend = await CheckBackendAsync(settings, cancellationToken);
        if (!backend.Ok) throw new InvalidOperationException(backend.Detail);

        var suffix = Path.GetExtension(sourceWindows).ToLowerInvariant();
        var faiWindows = sourceWindows;
        if (suffix is ".fa" or ".fasta" or ".fna")
        {
            var fastaWsl = PathBridge.WindowsToWsl(sourceWindows);
            var faidx = await RunWslAsync(
                settings.WslDistribution,
                RuntimeToolInvocation(settings, "samtools", "faidx", fastaWsl),
                cancellationToken);
            if (faidx.ExitCode != 0)
                throw new InvalidOperationException("FASTA konnte nicht mit samtools faidx indexiert werden.\n" + faidx.StdErr);
            faiWindows = sourceWindows + ".fai";
        }
        else if (suffix != ".fai")
        {
            throw new InvalidOperationException("Bitte eine FASTA (.fa/.fasta/.fna) oder deren .fai-Index auswählen.");
        }

        if (!File.Exists(faiWindows))
            throw new InvalidOperationException("Der erwartete FAI-Index wurde nicht gefunden: " + faiWindows);

        var faiSha256 = await Sha256FileAsync(faiWindows, cancellationToken);
        var referenceId = ReferenceIdFor(genomeBuild, faiSha256);

        var referenceDirWindows = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
            "ONTSeq", "references");
        Directory.CreateDirectory(referenceDirWindows);
        var lockWindows = Path.Combine(
            referenceDirWindows, ReferenceLockFileNameFor(genomeBuild, faiSha256));
        var temporaryLockWindows =
            lockWindows + "." + Guid.NewGuid().ToString("N") + ".tmp.json";
        var faiWsl = PathBridge.WindowsToWsl(faiWindows);
        var lockWsl = PathBridge.WindowsToWsl(lockWindows);
        var temporaryLockWsl = PathBridge.WindowsToWsl(temporaryLockWindows);

        try
        {
            var create = await RunWslAsync(
                settings.WslDistribution,
                BackendInvocation(settings,
                    "reference-lock",
                    "--fai", faiWsl,
                    "--reference-id", referenceId,
                    "--genome-build", genomeBuild,
                    "--require-canonical-assembly",
                    "--output", temporaryLockWsl),
                cancellationToken);
            if (create.ExitCode != 0)
                throw new InvalidOperationException(
                    "Reference-Lock konnte nicht erzeugt werden.\n" + create.StdErr);
            if (!File.Exists(temporaryLockWindows))
                throw new InvalidDataException(
                    "Reference-Lock wurde nicht am erwarteten temporären Ort geschrieben: " +
                    temporaryLockWindows);

            var validation = await CheckReferenceAsync(
                settings, temporaryLockWsl, genomeBuild, cancellationToken);
            if (!validation.Ok)
                throw new InvalidDataException(
                    "Der erzeugte Reference-Lock hat die Vollständigkeitsprüfung nicht " +
                    "bestanden.\n" + validation.Detail);

            PublishReferenceLockFile(temporaryLockWindows, lockWindows, faiSha256);
        }
        finally
        {
            if (File.Exists(temporaryLockWindows)) File.Delete(temporaryLockWindows);
        }

        settings.ReferenceLocksWsl[genomeBuild] = lockWsl;
        settings.SaveUserSettings();
        return lockWindows;
    }

    public static string ReferenceIdFor(string genomeBuild, string faiSha256) =>
        $"{genomeBuild}_LOCAL_{Sha256Prefix(faiSha256)}";

    public static string ReferenceLockFileNameFor(string genomeBuild, string faiSha256) =>
        $"{genomeBuild}.{Sha256Prefix(faiSha256)}.reference-lock.json";

    public static void PublishReferenceLockFile(
        string temporaryPath,
        string finalPath,
        string expectedFaiSha256)
    {
        using var document = JsonDocument.Parse(File.ReadAllText(temporaryPath));
        var recordedSha256 = document.RootElement.GetProperty("source_fai_sha256").GetString();
        if (!string.Equals(recordedSha256, expectedFaiSha256, StringComparison.Ordinal))
            throw new InvalidDataException(
                "Der erzeugte Reference-Lock enthält nicht den Fingerabdruck des gewählten FAI.");
        File.Move(temporaryPath, finalPath, overwrite: true);
    }

    private static string Sha256Prefix(string value)
    {
        if (value.Length != 64 || value.Any(character => !Uri.IsHexDigit(character)))
            throw new ArgumentException("Ein SHA256-Fingerabdruck mit 64 Hex-Zeichen wird erwartet.", nameof(value));
        return value[..16].ToLowerInvariant();
    }

    private static async Task<string> Sha256FileAsync(
        string path,
        CancellationToken cancellationToken)
    {
        await using var stream = new FileStream(
            path,
            FileMode.Open,
            FileAccess.Read,
            FileShare.Read,
            1024 * 1024,
            useAsync: true);
        return Convert.ToHexString(
            await SHA256.HashDataAsync(stream, cancellationToken)).ToLowerInvariant();
    }

    public async Task<string> RunSelfTestAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken)
    {
        var backend = await CheckBackendAsync(settings, cancellationToken);
        if (!backend.Ok) throw new InvalidOperationException(backend.Detail);
        if (string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
        {
            throw new InvalidOperationException(
                "Der vollständige System-Selbsttest benötigt die gebündelte ONTSeq-Runtime. " +
                "Bitte zuerst 'Runtime installieren' ausführen.");
        }

        var root = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
            "ONTSeq", "self-test", DateTime.Now.ToString("yyyyMMdd_HHmmss_fff"));
        Directory.CreateDirectory(root);
        var rootWsl = PathBridge.WindowsToWsl(root);
        var args = new List<string> { "system-smoke", "--output-dir", rootWsl };
        AddBundledPolicies(settings, args, includeCnv: true, includeCore034: false);
        var result = await RunWslAsync(
            settings.WslDistribution,
            BackendInvocation(settings, args.ToArray()),
            cancellationToken);
        File.WriteAllText(Path.Combine(root, "self-test.log.txt"), result.StdOut + Environment.NewLine + result.StdErr);
        if (result.ExitCode != 0)
            throw new InvalidOperationException(
                "ONTSeq Selbsttest ist fehlgeschlagen.\n" + result.StdErr);
        return root;
    }

    public void Start(
        DesktopSettings settings,
        string allowedRootWindows,
        string referenceLockWsl)
    {
        if (_process is { HasExited: false }) return;

        Directory.CreateDirectory(settings.OutputDirectoryWindows);
        StartProcess(settings, LegacyServiceArguments(settings, allowedRootWindows, referenceLockWsl).ToList());
    }

    internal static IReadOnlyList<string> LegacyServiceArguments(
        DesktopSettings settings,
        string allowedRootWindows,
        string referenceLockWsl)
    {
        var rootWsl = PathBridge.WindowsToWsl(allowedRootWindows);
        var outputWsl = PathBridge.WindowsToWsl(settings.OutputDirectoryWindows);
        var serviceArgs = new List<string>
        {
            "serve",
            "--reference-lock", referenceLockWsl,
            "--allow-root", rootWsl,
            "--output-dir", outputWsl,
            "--port", settings.Port.ToString(),
            "--no-browser"
        };
        if (!string.IsNullOrWhiteSpace(settings.AdaptiveTargetBedWsl))
        {
            serviceArgs.Add("--allow-root");
            serviceArgs.Add(PosixDirectoryName(settings.AdaptiveTargetBedWsl));
        }
        AddBundledPolicies(settings, serviceArgs, includeCnv: true, includeCore034: true);
        AddModkitArguments(settings, serviceArgs);
        return serviceArgs;
    }

    private void StartProcess(DesktopSettings settings, List<string> serviceArgs)
    {
        var args = BackendInvocation(settings, serviceArgs.ToArray());

        var psi = WslProcessStartInfo(settings.WslDistribution, args);
        _process = new Process { StartInfo = psi, EnableRaisingEvents = true };
        if (!_process.Start()) throw new InvalidOperationException("ONTSeq Backend konnte nicht gestartet werden.");
        _stdoutPump = CaptureLogAsync(_process.StandardOutput.BaseStream, _stdout);
        _stderrPump = CaptureLogAsync(_process.StandardError.BaseStream, _stderr);
    }

    public void StartProfile(
        DesktopSettings settings,
        string allowedRootWindows,
        string profileId) =>
        StartProfile(
            settings, allowedRootWindows, profileId, settings.Port, CreateServiceInstanceId());

    public void StartProfile(
        DesktopSettings settings,
        string allowedRootWindows,
        string profileId,
        int port) =>
        StartProfile(settings, allowedRootWindows, profileId, port, CreateServiceInstanceId());

    public void StartProfile(
        DesktopSettings settings,
        string allowedRootWindows,
        string profileId,
        int port,
        string instanceId)
    {
        _ = DesktopProfiles.Require(profileId);
        if (port is < 1 or > 65535)
            throw new ArgumentOutOfRangeException(nameof(port), port, "Ungültiger lokaler Dienstport.");
        if (_process is { HasExited: false }) return;

        settings.ApplyProfileDefaults();
        Directory.CreateDirectory(settings.OutputDirectoryWindows);
        StartProcess(
            settings, ProfileServiceArguments(settings, allowedRootWindows, port, instanceId).ToList());
    }

    internal static IReadOnlyList<string> ProfileServiceArguments(
        DesktopSettings settings,
        string allowedRootWindows,
        int port,
        string instanceId)
    {
        if (port is < 1 or > 65535)
            throw new ArgumentOutOfRangeException(nameof(port), port, "Ungültiger lokaler Dienstport.");
        if (!Guid.TryParseExact(instanceId, "N", out _) ||
            !string.Equals(instanceId, instanceId.ToLowerInvariant(), StringComparison.Ordinal))
            throw new ArgumentException(
                "Die Dienstinstanz muss als 32-stellige, kleingeschriebene GUID ohne Trennzeichen angegeben werden.",
                nameof(instanceId));
        var rootWsl = PathBridge.WindowsToWsl(allowedRootWindows);
        var outputWsl = PathBridge.WindowsToWsl(settings.OutputDirectoryWindows);
        var serviceArgs = new List<string>
        {
            "serve",
            "--resource-root", settings.ResourceRootWsl,
            "--allow-root", rootWsl,
            "--output-dir", outputWsl,
            "--port", port.ToString(),
            "--instance-id", instanceId,
            "--no-browser"
        };
        AddBundledPolicies(settings, serviceArgs, includeCnv: true, includeCore034: true);
        AddModkitArguments(settings, serviceArgs);
        return serviceArgs;
    }

    public bool HasExited => _process is null || _process.HasExited;

    public async ValueTask DisposeAsync()
    {
        if (_process is null) return;
        try
        {
            if (!_process.HasExited)
            {
                _process.Kill(entireProcessTree: true);
                await _process.WaitForExitAsync().WaitAsync(TimeSpan.FromSeconds(5));
            }
            await Task.WhenAll(_stdoutPump ?? Task.CompletedTask, _stderrPump ?? Task.CompletedTask)
                .WaitAsync(TimeSpan.FromSeconds(5));
        }
        catch
        {
        }
        finally
        {
            _process.Dispose();
            _process = null;
            _stdoutPump = null;
            _stderrPump = null;
        }
    }

    public static async Task<(int ExitCode, string StdOut, string StdErr)> RunWslAsync(
        string distribution,
        IReadOnlyList<string> command,
        CancellationToken cancellationToken)
    {
        var psi = WslProcessStartInfo(distribution, command);
        using var process = Process.Start(psi) ?? throw new InvalidOperationException("wsl.exe konnte nicht gestartet werden.");
        // Read bytes rather than the platform-default StreamReader: Linux commands emit
        // UTF-8, while wsl.exe's own bootstrap errors may be BOM-less UTF-16LE.
        var stdout = ReadWslOutputAsync(process.StandardOutput.BaseStream, cancellationToken);
        var stderr = ReadWslOutputAsync(process.StandardError.BaseStream, cancellationToken);
        await process.WaitForExitAsync(cancellationToken);
        return (process.ExitCode, await stdout, await stderr);
    }

    internal static ProcessStartInfo WslProcessStartInfo(string distribution, IReadOnlyList<string> command)
    {
        var psi = new ProcessStartInfo
        {
            FileName = "wsl.exe",
            UseShellExecute = false,
            CreateNoWindow = true,
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            StandardOutputEncoding = new UTF8Encoding(false, true),
            StandardErrorEncoding = new UTF8Encoding(false, true)
        };
        psi.ArgumentList.Add("-d");
        psi.ArgumentList.Add(distribution);
        psi.ArgumentList.Add("--");
        foreach (var item in command) psi.ArgumentList.Add(item);

        return psi;
    }

    private async Task CaptureLogAsync(Stream stream, StringBuilder destination)
    {
        try
        {
            await PumpWslOutputAsync(stream, text =>
            {
                lock (_outputLock) destination.Append(text);
            }, CancellationToken.None);
        }
        catch (DecoderFallbackException)
        {
            lock (_outputLock) destination.AppendLine("\nWSL-Ausgabe hat eine ungültige Textcodierung.");
        }
        catch (IOException) { /* A closed process pipe cannot supply further output. */ }
        catch (ObjectDisposedException) { /* The launcher has been disposed. */ }
    }

    internal static async Task<string> ReadWslOutputAsync(Stream stream, CancellationToken cancellationToken)
    {
        var output = new StringBuilder();
        await PumpWslOutputAsync(stream, text => output.Append(text), cancellationToken);
        return output.ToString();
    }

    private static async Task PumpWslOutputAsync(
        Stream stream, Action<string> append, CancellationToken cancellationToken)
    {
        var buffer = new byte[4096];
        var prefixLength = 0;
        while (prefixLength < 4)
        {
            var count = await stream.ReadAsync(buffer.AsMemory(prefixLength, 4 - prefixLength), cancellationToken);
            if (count == 0) break;
            prefixLength += count;
        }
        if (prefixLength == 0) return;

        Encoding encoding = new UTF8Encoding(false, true);
        var skip = 0;
        if (prefixLength >= 2 && buffer[0] == 0xff && buffer[1] == 0xfe)
        {
            encoding = new UnicodeEncoding(false, false, true);
            skip = 2;
        }
        else if (prefixLength >= 2 && buffer[0] == 0xfe && buffer[1] == 0xff)
        {
            encoding = new UnicodeEncoding(true, false, true);
            skip = 2;
        }
        else if (prefixLength >= 4 && buffer[0] != 0 && buffer[1] == 0 && buffer[2] != 0 && buffer[3] == 0)
            encoding = new UnicodeEncoding(false, false, true);
        else if (prefixLength >= 4 && buffer[0] == 0 && buffer[1] != 0 && buffer[2] == 0 && buffer[3] != 0)
            encoding = new UnicodeEncoding(true, false, true);
        else if (prefixLength >= 3 && buffer[0] == 0xef && buffer[1] == 0xbb && buffer[2] == 0xbf)
            skip = 3;

        // A decoder retains partial multibyte characters across pipe reads. No string
        // replacement or locale/code-page guessing is applied to already decoded text.
        var decoder = encoding.GetDecoder();
        var characters = new char[encoding.GetMaxCharCount(buffer.Length)];
        var characterCount = decoder.GetChars(buffer, skip, prefixLength - skip, characters, 0, false);
        if (characterCount > 0) append(new string(characters, 0, characterCount));
        while (true)
        {
            var count = await stream.ReadAsync(buffer, cancellationToken);
            characterCount = decoder.GetChars(buffer, 0, count, characters, 0, count == 0);
            if (characterCount > 0) append(new string(characters, 0, characterCount));
            if (count == 0) return;
        }
    }

    private static void AddBundledPolicies(
        DesktopSettings settings,
        List<string> args,
        bool includeCnv,
        bool includeCore034)
    {
        args.AddRange(BundledPolicyArguments(settings, includeCnv, includeCore034));
    }

    private static void AddModkitArguments(DesktopSettings settings, List<string> args)
    {
        var executable = DesktopSettings.ValidateModkitExecutableWsl(settings.ModkitExecutableWsl);
        if (executable is null) return;
        args.Add("--modkit");
        args.Add(executable);
    }

    internal static IReadOnlyList<string> BundledPolicyArguments(
        DesktopSettings settings,
        bool includeCnv,
        bool includeCore034)
    {
        var args = new List<string>();
        if (string.IsNullOrWhiteSpace(settings.RuntimeBinWsl)) return args;
        args.Add("--qc-policy");
        args.Add(RuntimeResource(settings, RuntimeAssets.QcPolicy));
        args.Add("--sniffles-policy");
        args.Add(RuntimeResource(settings, RuntimeAssets.SnifflesPolicy));
        if (includeCore034)
        {
            // These arguments belong to the runtime service parser. The engineering
            // system-smoke command deliberately has a smaller CLI contract.
            args.Add("--cutesv-policy");
            args.Add(RuntimeResource(settings, RuntimeAssets.CuteSvPolicy));
            args.Add("--sv-consensus-policy");
            args.Add(RuntimeResource(settings, RuntimeAssets.SvConsensusPolicy));
            args.Add("--sv-evidence-policy");
            args.Add(RuntimeResource(settings, RuntimeAssets.SvEvidencePolicy));
            args.Add("--target-coverage-policy");
            args.Add(RuntimeResource(settings, RuntimeAssets.TargetCoveragePolicy));
            args.Add("--components");
            args.Add(RuntimeResource(settings, RuntimeAssets.Components));
        }
        if (includeCnv)
        {
            args.Add("--cnv-policy");
            args.Add(RuntimeResource(settings, RuntimeAssets.CnvPolicy));
            args.Add("--qdnaseq-rscript");
            args.Add(RuntimeTool(settings, "Rscript"));
            args.Add("--qdnaseq-script");
            args.Add(RuntimeResource(settings, RuntimeAssets.QdnaSeqScript));
        }
        return args;
    }

    internal static IReadOnlyList<string> RequiredRuntimeFiles(DesktopSettings settings) =>
        RuntimeAssets.RequiredFiles.Select(asset => RuntimeResource(settings, asset)).ToArray();

    internal static IReadOnlyList<string> RequiredRuntimeTools(DesktopSettings settings) =>
        RuntimeAssets.RequiredTools.Select(tool => RuntimeTool(settings, tool)).ToArray();

    private static string BundledRuntimePrerequisiteCommand(DesktopSettings settings)
    {
        var checks = RequiredRuntimeFiles(settings)
            .Select(path => $"test -f {ShellQuote(path)}")
            .Concat(RequiredRuntimeTools(settings).Select(path => $"test -x {ShellQuote(path)}"));
        return string.Join(" && ", checks);
    }

    private static string RuntimeResource(DesktopSettings settings, string relative)
    {
        if (string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
            throw new InvalidOperationException("Gebündelte Runtime ist nicht konfiguriert.");
        var bin = settings.RuntimeBinWsl.TrimEnd('/');
        var root = bin.EndsWith("/bin", StringComparison.Ordinal) ? bin[..^4] : bin;
        return root + "/share/ontseq/" + relative.TrimStart('/');
    }

    private static string RuntimeTool(DesktopSettings settings, string tool)
    {
        if (string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
            throw new InvalidOperationException("Gebündelte Runtime ist nicht konfiguriert.");
        return settings.RuntimeBinWsl.TrimEnd('/') + "/" + tool;
    }

    private static string PosixDirectoryName(string path)
    {
        var normalized = path.Trim().TrimEnd('/');
        var separator = normalized.LastIndexOf('/');
        if (separator <= 0)
            throw new InvalidOperationException($"Kein absoluter WSL-Pfad: {path}");
        return normalized[..separator];
    }

    private static string ShellPathExpression(string path)
    {
        var normalized = DesktopSettings.NormalizeResourceRootWsl(path);
        if (!normalized.StartsWith("~/", StringComparison.Ordinal)) return ShellQuote(normalized);

        // Expand only the trusted $HOME prefix.  Escape every shell-significant character
        // in the user-configurable suffix so it remains one path inside the home directory.
        var suffix = normalized[2..]
            .Replace("\\", "\\\\", StringComparison.Ordinal)
            .Replace("\"", "\\\"", StringComparison.Ordinal)
            .Replace("$", "\\$", StringComparison.Ordinal)
            .Replace("`", "\\`", StringComparison.Ordinal);
        return $"\"$HOME/{suffix}\"";
    }

    internal static IReadOnlyList<string> BackendInvocation(DesktopSettings settings, params string[] args)
    {
        var command = new List<string> { "env", "PYTHONIOENCODING=utf-8", "PYTHONUTF8=1" };
        if (!string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
        {
            command.Add($"PATH={settings.RuntimeBinWsl}:{BaseLinuxPath}");
        }
        command.Add(settings.BackendCommand);
        command.AddRange(args);
        return command;
    }

    private static IReadOnlyList<string> RuntimeToolInvocation(
        DesktopSettings settings,
        string tool,
        params string[] args)
    {
        var executable = string.IsNullOrWhiteSpace(settings.RuntimeBinWsl)
            ? tool
            : RuntimeTool(settings, tool);
        var command = new List<string>();
        if (!string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
        {
            command.Add("env");
            command.Add($"PATH={settings.RuntimeBinWsl}:{BaseLinuxPath}");
        }
        command.Add(executable);
        command.AddRange(args);
        return command;
    }

    private static string ShellQuote(string value) => "'" + value.Replace("'", "'\"'\"'") + "'";
}

public static class PathBridge
{
    public static string WindowsToWsl(string path)
    {
        var selected = path.Trim().Trim('"');
        if (selected.StartsWith("\\\\", StringComparison.Ordinal) || selected.StartsWith("//", StringComparison.Ordinal))
            throw new InvalidOperationException(
                "UNC- und Gerätepfade werden nicht automatisch nach WSL übersetzt. Wähle einen vollständigen Laufwerkspfad.");
        // Reject relative and Linux paths before GetFullPath can silently reinterpret
        // them against the Desktop's current Windows directory or current drive.
        if (selected.Any(char.IsControl) || selected.Length < 3 ||
            !char.IsAsciiLetter(selected[0]) || selected[1] != ':' ||
            (selected[2] != '\\' && selected[2] != '/'))
            throw new InvalidOperationException("Erwartet wird ein vollständiger Windows-Laufwerkspfad, z. B. E:\\Daten.");
        var full = Path.GetFullPath(selected);
        var drive = char.ToLowerInvariant(full[0]);
        var remainder = full[3..].Replace('\\', '/');
        // This is the default WSL drive mapping, not evidence that the drive is mounted
        // or accessible. The named distribution's prerequisite checks establish that.
        return string.IsNullOrEmpty(remainder) ? $"/mnt/{drive}" : $"/mnt/{drive}/{remainder}";
    }
}
