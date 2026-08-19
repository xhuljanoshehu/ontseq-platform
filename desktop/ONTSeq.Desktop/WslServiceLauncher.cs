using System.Diagnostics;
using System.Text;

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
        CancellationToken cancellationToken)
    {
        var status = await CheckWslAsync(settings, cancellationToken);
        if (!status.Ok) throw new InvalidOperationException(status.Detail);

        var backend = await CheckBackendAsync(settings, cancellationToken);
        if (!backend.Ok) throw new InvalidOperationException(backend.Detail);

        var rootWsl = PathBridge.WindowsToWsl(allowedRootWindows);
        var outputWsl = PathBridge.WindowsToWsl(settings.OutputDirectoryWindows);
        var checks = $"test -d {ShellQuote(rootWsl)} && test -f {ShellQuote(referenceLockWsl)} && mkdir -p {ShellQuote(outputWsl)}";
        var check = await RunWslAsync(settings.WslDistribution, ["sh", "-lc", checks], cancellationToken);
        if (check.ExitCode != 0)
        {
            throw new InvalidOperationException(
                "Der Windows-Speicher oder das Referenzpaket ist in WSL nicht erreichbar. " +
                "Bei einem Netzlaufwerk wie P: muss dieses in WSL als drvfs eingebunden sein.\n" +
                $"Erwarteter WSL-Pfad: {rootWsl}\nReferenz: {referenceLockWsl}\n{check.StdErr}");
        }
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
            return result.ExitCode == 0
                ? (true, $"ONTSeq Backend gefunden: {settings.BackendCommand}")
                : (false, "ONTSeq Backend ist in WSL nicht einsatzbereit. " + result.StdErr);
        }
        catch (Exception error)
        {
            return (false, "ONTSeq Backend fehlt oder kann nicht gestartet werden: " + error.Message);
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

        var wsl = await CheckWslAsync(settings, cancellationToken);
        if (!wsl.Ok) throw new InvalidOperationException(wsl.Detail);

        var homeResult = await RunWslAsync(
            settings.WslDistribution, ["sh", "-lc", "printf %s \"$HOME\""], cancellationToken);
        if (homeResult.ExitCode != 0 || string.IsNullOrWhiteSpace(homeResult.StdOut))
            throw new InvalidOperationException("WSL-Home-Verzeichnis konnte nicht bestimmt werden. " + homeResult.StdErr);

        var home = homeResult.StdOut.Trim();
        var target = home + "/.local/share/ontseq/runtime-v0.1.1";
        var bin = target + "/bin";
        var archiveWsl = PathBridge.WindowsToWsl(runtimeArchiveWindows);
        var command =
            $"rm -rf {ShellQuote(target)} && mkdir -p {ShellQuote(target)} && " +
            $"tar -xzf {ShellQuote(archiveWsl)} -C {ShellQuote(target)} && " +
            $"{ShellQuote(bin + "/conda-unpack")} && " +
            $"env PATH={ShellQuote(bin + ":" + BaseLinuxPath)} {ShellQuote(bin + "/ontseq")} --help >/dev/null";
        var install = await RunWslAsync(
            settings.WslDistribution, ["sh", "-lc", command], cancellationToken);
        if (install.ExitCode != 0)
            throw new InvalidOperationException("ONTSeq Linux-Runtime konnte nicht installiert werden.\n" + install.StdErr);

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

        var referenceDirWindows = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
            "ONTSeq", "references");
        Directory.CreateDirectory(referenceDirWindows);
        var lockWindows = Path.Combine(referenceDirWindows, $"{genomeBuild}.reference-lock.json");
        var faiWsl = PathBridge.WindowsToWsl(faiWindows);
        var lockWsl = PathBridge.WindowsToWsl(lockWindows);
        var referenceId = $"{genomeBuild}_LOCAL_{DateTime.UtcNow:yyyyMMdd}";

        var create = await RunWslAsync(
            settings.WslDistribution,
            BackendInvocation(settings,
                "reference-lock",
                "--fai", faiWsl,
                "--reference-id", referenceId,
                "--genome-build", genomeBuild,
                "--output", lockWsl),
            cancellationToken);
        if (create.ExitCode != 0)
            throw new InvalidOperationException("Reference-Lock konnte nicht erzeugt werden.\n" + create.StdErr);
        if (!File.Exists(lockWindows))
            throw new InvalidOperationException("Reference-Lock wurde nicht am erwarteten Ort geschrieben: " + lockWindows);

        settings.ReferenceLocksWsl[genomeBuild] = lockWsl;
        settings.SaveUserSettings();
        return lockWindows;
    }

    public async Task<string> RunSelfTestAsync(
        DesktopSettings settings,
        CancellationToken cancellationToken)
    {
        var backend = await CheckBackendAsync(settings, cancellationToken);
        if (!backend.Ok) throw new InvalidOperationException(backend.Detail);

        var root = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
            "ONTSeq", "self-test", DateTime.Now.ToString("yyyyMMdd_HHmmss"));
        Directory.CreateDirectory(root);
        var rootWsl = PathBridge.WindowsToWsl(root);
        var result = await RunWslAsync(
            settings.WslDistribution,
            BackendInvocation(settings, "local-smoke", "--output-dir", rootWsl),
            cancellationToken);
        File.WriteAllText(Path.Combine(root, "self-test.log.txt"), result.StdOut + Environment.NewLine + result.StdErr);
        if (result.ExitCode != 0)
            throw new InvalidOperationException("ONTSeq Selbsttest ist fehlgeschlagen.\n" + result.StdErr);
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

        var args = BackendInvocation(
            settings,
            "serve",
            "--reference-lock", referenceLockWsl,
            "--allow-root", rootWsl,
            "--output-dir", outputWsl,
            "--port", settings.Port.ToString(),
            "--no-browser");

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
            // App shutdown must not be blocked by a backend that is already disappearing.
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
