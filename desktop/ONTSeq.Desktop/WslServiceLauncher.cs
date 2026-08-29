using System.Diagnostics;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;

namespace ONTSeq.Desktop;

public sealed class WslServiceLauncher : IAsyncDisposable
{
    private const string BaseLinuxPath = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin";
    private const string ReleaseVersion = "0.5.3";

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
    public const string HematologyKnowledgeBundleId = "HEMATOLOGY_v3";
    public const string AmlAdaptivePanelBundleId = "AML_AS_111_GRCh38_v1";
    public static IReadOnlyList<string> ManagedGrch38ResourceBundleIds { get; } =
        Array.AsReadOnly(new[] {
            Grch38ReferenceBundleId,
            HematologyKnowledgeBundleId,
            AmlAdaptivePanelBundleId
        });
    private Process? _process;
    private readonly StringBuilder _stderr = new();
    private readonly StringBuilder _stdout = new();

    public string DiagnosticLog => $"{_stdout}\n{_stderr}".Trim();

    public async Task VerifyPrerequisitesAsync(
        DesktopSettings settings,
        string allowedRootWindows,
        string referenceLockWsl,
        string expectedGenomeBuild,
        CancellationToken cancellationToken)
    {
        var status = await CheckWslAsync(settings, cancellationToken);
        if (!status.Ok) throw new InvalidOperationException(status.Detail);

        var backend = await CheckBackendAsync(settings, cancellationToken);
        if (!backend.Ok) throw new InvalidOperationException(backend.Detail);

        var rootWsl = PathBridge.WindowsToWsl(allowedRootWindows);
        var outputWsl = PathBridge.WindowsToWsl(settings.OutputDirectoryWindows);
        var checks = $"test -d {ShellQuote(rootWsl)} && test -f {ShellQuote(referenceLockWsl)} && mkdir -p {ShellQuote(outputWsl)}";
        if (!string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
        {
            checks += " && " + BundledRuntimePrerequisiteCommand(settings);
        }
        if (!string.IsNullOrWhiteSpace(settings.AdaptiveTargetBedWsl))
            checks += $" && test -s {ShellQuote(settings.AdaptiveTargetBedWsl)}";

        var check = await RunWslAsync(settings.WslDistribution, ["sh", "-lc", checks], cancellationToken);
        if (check.ExitCode != 0)
        {
            throw new InvalidOperationException(
                "Der Windows-Speicher, die Runtime-Ressourcen oder das Referenzpaket sind in WSL nicht erreichbar. " +
                "Bei einem Netzlaufwerk wie P: muss dieses in WSL als drvfs eingebunden sein.\n" +
                $"Erwarteter WSL-Pfad: {rootWsl}\nReferenz: {referenceLockWsl}\n{check.StdErr}");
        }

        var reference = await CheckReferenceAsync(
            settings, referenceLockWsl, expectedGenomeBuild, cancellationToken);
        if (!reference.Ok)
            throw new InvalidDataException(
                $"Die gespeicherte {expectedGenomeBuild}-Referenz ist ungültig oder " +
                "unvollständig. Öffne 'System einrichten' und wähle den vollständigen " +
                $"FAI-Index der BAM-Referenz.\n{reference.Detail}");
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
                !serviceCapability.StdOut.Contains("--resource-root", StringComparison.Ordinal))
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
                    "GRCh38-Bundles. Bitte 'Runtime installieren' erneut ausführen.");
            }

            return (true, $"ONTSeq Backend v{ReleaseVersion} gefunden: {settings.BackendCommand}");
        }
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

    public async Task<(bool Ok, string Detail)> CheckResourceBundlesAsync(
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
        return InterpretResourceStatus(result.ExitCode, result.StdOut, result.StdErr);
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
        string stderr)
    {
        if (exitCode != 0)
            return (false, (string.IsNullOrWhiteSpace(stderr) ? stdout : stderr).Trim());
        try
        {
            using var document = JsonDocument.Parse(stdout);
            var root = document.RootElement;
            var referenceReady = root.TryGetProperty("references", out var references) &&
                references.ValueKind == JsonValueKind.Array &&
                references.EnumerateArray().Any(item =>
                    item.TryGetProperty("bundle_id", out var id) &&
                    string.Equals(id.GetString(), Grch38ReferenceBundleId, StringComparison.Ordinal) &&
                    item.TryGetProperty("valid", out var valid) && valid.GetBoolean());
            var profiles = new HashSet<string>(StringComparer.Ordinal);
            if (root.TryGetProperty("profiles", out var profileArray) &&
                profileArray.ValueKind == JsonValueKind.Array)
            {
                foreach (var item in profileArray.EnumerateArray())
                {
                    if (item.ValueKind == JsonValueKind.String && item.GetString() is { } profileId)
                        profiles.Add(profileId);
                }
            }
            var requiredProfilesReady = DesktopProfiles.Supported.All(
                profile => profiles.Contains(profile.ProfileId));
            if (referenceReady && requiredProfilesReady)
            {
                return (
                    true,
                    "GRCh38-Ressourcenfamilie gültig (" +
                    string.Join(", ", ManagedGrch38ResourceBundleIds) + ") · Profile: " +
                    string.Join(", ", DesktopProfiles.Supported.Select(item => item.ProfileId)));
            }

            var missing = new List<string>();
            if (!referenceReady) missing.Add(Grch38ReferenceBundleId);
            missing.AddRange(DesktopProfiles.Supported
                .Where(profile => !profiles.Contains(profile.ProfileId))
                .Select(profile => profile.ProfileId));
            return (false, "Nicht einsatzbereit: " + string.Join(", ", missing));
        }
        catch (JsonException error)
        {
            return (false, "Ungültige JSON-Antwort von 'ontseq references status': " + error.Message);
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
        _ = DesktopProfiles.Require(profileId);
        settings.ApplyProfileDefaults();

        var status = await CheckWslAsync(settings, cancellationToken);
        if (!status.Ok) throw new InvalidOperationException(status.Detail);
        var backend = await CheckBackendAsync(settings, cancellationToken);
        if (!backend.Ok) throw new InvalidOperationException(backend.Detail);

        var rootWsl = PathBridge.WindowsToWsl(allowedRootWindows);
        var outputWsl = PathBridge.WindowsToWsl(settings.OutputDirectoryWindows);
        var checks = new List<string>
        {
            $"test -d {ShellQuote(rootWsl)}",
            $"mkdir -p {ShellQuote(outputWsl)}",
            $"test -d {ShellPathExpression(settings.ResourceRootWsl)}"
        };
        if (!string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
            checks.Add(BundledRuntimePrerequisiteCommand(settings));
        var check = await RunWslAsync(
            settings.WslDistribution,
            ["sh", "-lc", string.Join(" && ", checks)],
            cancellationToken);
        if (check.ExitCode != 0)
        {
            throw new InvalidOperationException(
                "BAM-Speicher, Ausgabeverzeichnis, Resource-Root oder Linux-Runtime sind in WSL " +
                "nicht vollständig erreichbar. " +
                "Öffne 'System einrichten', installiere das GRCh38-Bundle und prüfe bei " +
                "Netzlaufwerken die drvfs-Einbindung.\n" +
                $"BAM-Root: {rootWsl}\nResource-Root: {settings.ResourceRootWsl}\n{check.StdErr}");
        }

        var statusReport = await CheckResourceBundlesAsync(settings, cancellationToken);
        if (!statusReport.Ok)
            throw new InvalidOperationException(
                "Die manifestierten GRCh38-Profile sind nicht vollständig einsatzbereit. " +
                "Nutze in 'System einrichten' die vollständige Ressourcen-Reparatur für " +
                $"{Grch38ReferenceBundleId}, {HematologyKnowledgeBundleId} und " +
                $"{AmlAdaptivePanelBundleId}.\n" +
                statusReport.Detail);
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

        var wsl = await CheckWslAsync(settings, cancellationToken);
        if (!wsl.Ok) throw new InvalidOperationException(wsl.Detail);

        var homeResult = await RunWslAsync(
            settings.WslDistribution, ["sh", "-lc", "printf %s \"$HOME\""], cancellationToken);
        if (homeResult.ExitCode != 0 || string.IsNullOrWhiteSpace(homeResult.StdOut))
            throw new InvalidOperationException("WSL-Home-Verzeichnis konnte nicht bestimmt werden. " + homeResult.StdErr);

        var home = homeResult.StdOut.Trim();
        var target = home + $"/.local/share/ontseq/runtime-v{ReleaseVersion}";
        var bin = target + "/bin";
        var runtimePath = bin + ":" + BaseLinuxPath;
        var archiveWsl = PathBridge.WindowsToWsl(runtimeArchiveWindows);
        var installedRuntime = new DesktopSettings { RuntimeBinWsl = bin };
        var command =
            $"rm -rf {ShellQuote(target)} && mkdir -p {ShellQuote(target)} && " +
            $"tar -xzf {ShellQuote(archiveWsl)} -C {ShellQuote(target)} && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/conda-unpack")} && " +
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
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--resource-root'";
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
        StartProcess(settings, serviceArgs);
    }

    private void StartProcess(DesktopSettings settings, List<string> serviceArgs)
    {
        var args = BackendInvocation(settings, serviceArgs.ToArray());

        var psi = new ProcessStartInfo
        {
            FileName = "wsl.exe",
            UseShellExecute = false,
            CreateNoWindow = true,
            RedirectStandardOutput = true,
            RedirectStandardError = true
        };
        psi.ArgumentList.Add("-d");
        psi.ArgumentList.Add(settings.WslDistribution);
        psi.ArgumentList.Add("--");
        foreach (var item in args) psi.ArgumentList.Add(item);

        _process = new Process { StartInfo = psi, EnableRaisingEvents = true };
        _process.OutputDataReceived += (_, e) => { if (e.Data is not null) _stdout.AppendLine(e.Data); };
        _process.ErrorDataReceived += (_, e) => { if (e.Data is not null) _stderr.AppendLine(e.Data); };
        if (!_process.Start()) throw new InvalidOperationException("ONTSeq Backend konnte nicht gestartet werden.");
        _process.BeginOutputReadLine();
        _process.BeginErrorReadLine();
    }

    public void StartProfile(
        DesktopSettings settings,
        string allowedRootWindows,
        string profileId)
    {
        _ = DesktopProfiles.Require(profileId);
        if (_process is { HasExited: false }) return;

        settings.ApplyProfileDefaults();
        Directory.CreateDirectory(settings.OutputDirectoryWindows);
        var rootWsl = PathBridge.WindowsToWsl(allowedRootWindows);
        var outputWsl = PathBridge.WindowsToWsl(settings.OutputDirectoryWindows);
        var serviceArgs = new List<string>
        {
            "serve",
            "--resource-root", settings.ResourceRootWsl,
            "--allow-root", rootWsl,
            "--output-dir", outputWsl,
            "--port", settings.Port.ToString(),
            "--no-browser"
        };
        AddBundledPolicies(settings, serviceArgs, includeCnv: true, includeCore034: true);
        StartProcess(settings, serviceArgs);
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
        }
        catch
        {
        }
        finally
        {
            _process.Dispose();
            _process = null;
        }
    }

    public static async Task<(int ExitCode, string StdOut, string StdErr)> RunWslAsync(
        string distribution,
        IReadOnlyList<string> command,
        CancellationToken cancellationToken)
    {
        var psi = new ProcessStartInfo
        {
            FileName = "wsl.exe",
            UseShellExecute = false,
            CreateNoWindow = true,
            RedirectStandardOutput = true,
            RedirectStandardError = true
        };
        psi.ArgumentList.Add("-d");
        psi.ArgumentList.Add(distribution);
        psi.ArgumentList.Add("--");
        foreach (var item in command) psi.ArgumentList.Add(item);

        using var process = Process.Start(psi) ?? throw new InvalidOperationException("wsl.exe konnte nicht gestartet werden.");
        var stdout = process.StandardOutput.ReadToEndAsync(cancellationToken);
        var stderr = process.StandardError.ReadToEndAsync(cancellationToken);
        await process.WaitForExitAsync(cancellationToken);
        return (process.ExitCode, await stdout, await stderr);
    }

    private static void AddBundledPolicies(
        DesktopSettings settings,
        List<string> args,
        bool includeCnv,
        bool includeCore034)
    {
        args.AddRange(BundledPolicyArguments(settings, includeCnv, includeCore034));
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

    private static IReadOnlyList<string> BackendInvocation(DesktopSettings settings, params string[] args)
    {
        var command = new List<string>();
        if (!string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
        {
            command.Add("env");
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
        var full = Path.GetFullPath(path.Trim().Trim('"'));
        if (full.StartsWith("\\\\", StringComparison.Ordinal))
            throw new InvalidOperationException(
                "UNC-Pfade werden nicht automatisch nach WSL übersetzt. Das Netzlaufwerk muss zuerst als drvfs gemountet werden.");
        if (full.Length < 3 || full[1] != ':' || (full[2] != '\\' && full[2] != '/'))
            throw new InvalidOperationException($"Kein unterstützter Windows-Pfad: {full}");
        var drive = char.ToLowerInvariant(full[0]);
        var remainder = full[3..].Replace('\\', '/');
        return string.IsNullOrEmpty(remainder) ? $"/mnt/{drive}" : $"/mnt/{drive}/{remainder}";
    }
}
