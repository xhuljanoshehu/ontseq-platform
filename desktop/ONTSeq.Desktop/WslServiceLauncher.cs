using System.Diagnostics;
using System.Text;

namespace ONTSeq.Desktop;

public sealed class WslServiceLauncher : IAsyncDisposable
{
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
        var status = await RunWslAsync(settings.WslDistribution, ["sh", "-lc", "printf ready"], cancellationToken);
        if (status.ExitCode != 0 || !status.StdOut.Contains("ready", StringComparison.Ordinal))
            throw new InvalidOperationException(
                "WSL2 ist nicht einsatzbereit. Windows muss WSL2 und die konfigurierte Linux-Distribution bereitstellen.\n" + status.StdErr);

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

    public void Start(
        DesktopSettings settings,
        string allowedRootWindows,
        string referenceLockWsl)
    {
        if (_process is { HasExited: false }) return;

        Directory.CreateDirectory(settings.OutputDirectoryWindows);
        var rootWsl = PathBridge.WindowsToWsl(allowedRootWindows);
        var outputWsl = PathBridge.WindowsToWsl(settings.OutputDirectoryWindows);

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
        psi.ArgumentList.Add(settings.BackendCommand);
        psi.ArgumentList.Add("serve");
        psi.ArgumentList.Add("--reference-lock");
        psi.ArgumentList.Add(referenceLockWsl);
        psi.ArgumentList.Add("--allow-root");
        psi.ArgumentList.Add(rootWsl);
        psi.ArgumentList.Add("--output-dir");
        psi.ArgumentList.Add(outputWsl);
        psi.ArgumentList.Add("--port");
        psi.ArgumentList.Add(settings.Port.ToString());
        psi.ArgumentList.Add("--no-browser");

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

    private static async Task<(int ExitCode, string StdOut, string StdErr)> RunWslAsync(
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
