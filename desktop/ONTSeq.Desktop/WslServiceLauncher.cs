using System.Diagnostics;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;

namespace ONTSeq.Desktop;

public sealed class WslServiceLauncher : IAsyncDisposable
{
    private const string BaseLinuxPath = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin";
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
            checks += $" && test -f {ShellQuote(RuntimeResource(settings, "configs/qc/defaults.yaml"))}" +
                      $" && test -f {ShellQuote(RuntimeResource(settings, "configs/qc/adaptive_target_coverage.technical.yaml"))}" +
                      $" && test -f {ShellQuote(RuntimeResource(settings, "configs/components/default.yaml"))}" +
                      $" && test -f {ShellQuote(RuntimeResource(settings, "configs/sv/sniffles2.conservative.technical.yaml"))}" +
                      $" && test -f {ShellQuote(RuntimeResource(settings, "configs/cnv/qdnaseq_ace.technical.yaml"))}" +
                      $" && test -f {ShellQuote(RuntimeResource(settings, "scripts/run_qdnaseq_ace.R"))}";
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
                    "installieren' ausführen, um sie auf Desktop/Core v0.3.4 zu aktualisieren.");

            var serviceCapability = await RunWslAsync(
                settings.WslDistribution,
                BackendInvocation(settings, "serve", "--help"),
                cancellationToken);
            if (serviceCapability.ExitCode != 0 ||
                !serviceCapability.StdOut.Contains("--target-coverage-policy", StringComparison.Ordinal) ||
                !serviceCapability.StdOut.Contains("--components", StringComparison.Ordinal))
            {
                return (
                    false,
                    "Die installierte ONTSeq Runtime enthält nicht den vollständigen v0.3.4-" +
                    "Desktop-Vertrag für Target Coverage und Komponentenauswahl. Bitte " +
                    "'Runtime installieren' erneut ausführen.");
            }

            return (true, $"ONTSeq Backend v0.3.4 gefunden: {settings.BackendCommand}");
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
        var target = home + "/.local/share/ontseq/runtime-v0.3.4";
        var bin = target + "/bin";
        var runtimePath = bin + ":" + BaseLinuxPath;
        var archiveWsl = PathBridge.WindowsToWsl(runtimeArchiveWindows);
        var qc = target + "/share/ontseq/configs/qc/defaults.yaml";
        var targetCoverage = target + "/share/ontseq/configs/qc/adaptive_target_coverage.technical.yaml";
        var components = target + "/share/ontseq/configs/components/default.yaml";
        var sniffles = target + "/share/ontseq/configs/sv/sniffles2.conservative.technical.yaml";
        var cnvPolicy = target + "/share/ontseq/configs/cnv/qdnaseq_ace.technical.yaml";
        var cnvScript = target + "/share/ontseq/scripts/run_qdnaseq_ace.R";
        var command =
            $"rm -rf {ShellQuote(target)} && mkdir -p {ShellQuote(target)} && " +
            $"tar -xzf {ShellQuote(archiveWsl)} -C {ShellQuote(target)} && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/conda-unpack")} && " +
            $"test -f {ShellQuote(qc)} && test -f {ShellQuote(targetCoverage)} && " +
            $"test -f {ShellQuote(components)} && test -f {ShellQuote(sniffles)} && " +
            $"test -f {ShellQuote(cnvPolicy)} && test -f {ShellQuote(cnvScript)} && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/Rscript")} -e " +
            ShellQuote("stopifnot(requireNamespace('QDNAseq',quietly=TRUE), requireNamespace('QDNAseq.hg19',quietly=TRUE), requireNamespace('ACE',quietly=TRUE))") +
            " && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} validate-reference --help >/dev/null && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--target-coverage-policy' && " +
            $"env PATH={ShellQuote(runtimePath)} {ShellQuote(bin + "/ontseq")} serve --help | grep -q -- '--components'";
        var install = await RunWslAsync(
            settings.WslDistribution, ["sh", "-lc", command], cancellationToken);
        if (install.ExitCode != 0)
            throw new InvalidOperationException("ONTSeq Linux-Runtime v0.3.4 konnte nicht installiert werden.\n" + install.StdErr);

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
        if (string.IsNullOrWhiteSpace(settings.RuntimeBinWsl)) return;
        args.Add("--qc-policy");
        args.Add(RuntimeResource(settings, "configs/qc/defaults.yaml"));
        args.Add("--sniffles-policy");
        args.Add(RuntimeResource(settings, "configs/sv/sniffles2.conservative.technical.yaml"));
        if (includeCore034)
        {
            args.Add("--target-coverage-policy");
            args.Add(RuntimeResource(settings, "configs/qc/adaptive_target_coverage.technical.yaml"));
            args.Add("--components");
            args.Add(RuntimeResource(settings, "configs/components/default.yaml"));
        }
        if (!includeCnv) return;
        args.Add("--cnv-policy");
        args.Add(RuntimeResource(settings, "configs/cnv/qdnaseq_ace.technical.yaml"));
        args.Add("--qdnaseq-rscript");
        args.Add(settings.RuntimeBinWsl.TrimEnd('/') + "/Rscript");
        args.Add("--qdnaseq-script");
        args.Add(RuntimeResource(settings, "scripts/run_qdnaseq_ace.R"));
    }

    private static string RuntimeResource(DesktopSettings settings, string relative)
    {
        if (string.IsNullOrWhiteSpace(settings.RuntimeBinWsl))
            throw new InvalidOperationException("Gebündelte Runtime ist nicht konfiguriert.");
        var root = settings.RuntimeBinWsl.EndsWith("/bin", StringComparison.Ordinal)
            ? settings.RuntimeBinWsl[..^4]
            : settings.RuntimeBinWsl.TrimEnd('/');
        return root + "/share/ontseq/" + relative.TrimStart('/');
    }

    private static string PosixDirectoryName(string path)
    {
        var normalized = path.Trim().TrimEnd('/');
        var separator = normalized.LastIndexOf('/');
        if (separator <= 0)
            throw new InvalidOperationException($"Kein absoluter WSL-Pfad: {path}");
        return normalized[..separator];
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
            : settings.RuntimeBinWsl.TrimEnd('/') + "/" + tool;
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
