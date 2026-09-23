using System.Diagnostics;
using System.Text.Json;

namespace ONTSeq.Desktop;

/// <summary>Repairs a proven stale local data drive after the owned backend has exited.</summary>
public static class WslDriveRecovery
{
    public static async Task<IReadOnlyList<string>> CheckAsync(
        DesktopSettings settings, IEnumerable<string> paths, CancellationToken cancellationToken)
    {
        var repaired = new List<string>();
        var drives = paths.Where(path => !string.IsNullOrWhiteSpace(path)).Select(path =>
        {
            _ = PathBridge.WindowsToWsl(path); // Reject UNC, relative and malformed paths.
            return char.ToUpperInvariant(path[0]);
        }).Distinct().ToArray();
        using var stream = typeof(WslDriveRecovery).Assembly.GetManifestResourceStream(
            "ONTSeq.Desktop.wsl_drive_recovery.py") ?? throw new InvalidDataException("Laufwerksprüfung fehlt im Paket.");
        using var reader = new StreamReader(stream);
        var script = await reader.ReadToEndAsync(cancellationToken);
        var user = (await RunAsync(settings.WslDistribution, ["/usr/bin/id", "-u"], false, cancellationToken)).Trim();
        var group = (await RunAsync(settings.WslDistribution, ["/usr/bin/id", "-g"], false, cancellationToken)).Trim();
        if (!int.TryParse(user, out var uid) || !int.TryParse(group, out var gid) || uid <= 0 || gid <= 0)
            throw new InvalidOperationException("WSL muss einen normalen Standardbenutzer verwenden.");
        foreach (var drive in drives)
        {
            // C: and network drives are checked by ordinary prerequisites, never remounted here.
            if (drive < 'D' || drive > 'Z') continue;
            var windowsDrive = new DriveInfo(drive + @":\");
            if (!windowsDrive.IsReady)
                throw new IOException($"{drive}: ist bereits unter Windows nicht erreichbar. Laufwerk verbinden und erneut versuchen.");
            if (windowsDrive.DriveType is not (DriveType.Fixed or DriveType.Removable)) continue;
            var args = new[] { "/usr/bin/python3", "-I", "-c", script, drive.ToString(), user, group, "check" };
            using var status = JsonDocument.Parse(await RunAsync(settings.WslDistribution, args, false, cancellationToken));
            if (status.RootElement.GetProperty("state").GetString() != "stale") continue;
            args[^1] = "repair";
            using var result = JsonDocument.Parse(await RunAsync(settings.WslDistribution, args, true, cancellationToken));
            if (result.RootElement.GetProperty("state").GetString() is not ("repaired" or "healthy"))
                throw new IOException($"Die Reparatur von {drive}: wurde nicht bestätigt.");
            repaired.Add(drive + ":");
        }
        return repaired;
    }

    private static async Task<string> RunAsync(string distribution, IReadOnlyList<string> command,
        bool asRoot, CancellationToken cancellationToken)
    {
        var info = WslServiceLauncher.WslProcessStartInfo(distribution, command);
        if (asRoot)
        {
            info.ArgumentList.Insert(2, "-u");
            info.ArgumentList.Insert(3, "root");
        }
        using var timeout = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
        timeout.CancelAfter(TimeSpan.FromSeconds(40));
        using var process = Process.Start(info) ?? throw new IOException("WSL-Prüfung konnte nicht starten.");
        var stdout = WslServiceLauncher.ReadWslOutputAsync(process.StandardOutput.BaseStream, timeout.Token);
        var stderr = WslServiceLauncher.ReadWslOutputAsync(process.StandardError.BaseStream, timeout.Token);
        try
        {
            await process.WaitForExitAsync(timeout.Token);
            var output = await stdout;
            var error = await stderr;
            if (process.ExitCode != 0)
                throw new IOException("Laufwerksprüfung/Reparatur nicht abgeschlossen. " + output + error);
            return output;
        }
        finally
        {
            // Only this short-lived helper is owned here; no distribution/process-name kill.
            if (!process.HasExited) process.Kill(entireProcessTree: true);
            try { await Task.WhenAll(stdout, stderr); } catch (OperationCanceledException) { }
        }
    }
}
