using System.Security.Cryptography;
using System.Text.RegularExpressions;

namespace ONTSeq.Desktop;

public sealed record RuntimePackage(string ArchivePath, string WheelPath)
{
    public static async Task<RuntimePackage> VerifyAsync(
        string archivePath, string coreVersion, CancellationToken cancellationToken)
    {
        if (!Regex.IsMatch(coreVersion, @"\A[0-9]+\.[0-9]+\.[0-9]+\z"))
            throw new ArgumentException("Ungültige Core-Version.", nameof(coreVersion));
        var archive = Path.GetFullPath(archivePath);
        if (Path.GetFileName(archive) != "ontseq-linux-runtime.tar.gz")
            throw new InvalidDataException("Unerwarteter Dateiname für die gebündelte Runtime.");
        var directory = Path.GetDirectoryName(archive)!;
        var wheel = Path.Combine(directory, $"ontseq_platform-{coreVersion}-py3-none-any.whl");
        var checksumsPath = Path.Combine(directory, "SHA256SUMS");
        var checksums = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (var line in await File.ReadAllLinesAsync(checksumsPath, cancellationToken))
        {
            if (string.IsNullOrWhiteSpace(line) || line.StartsWith('#')) continue;
            var match = Regex.Match(line, @"\A(?<hash>[0-9a-fA-F]{64}) [ *](?<name>[^\r\n]+)\z");
            if (!match.Success)
                throw new InvalidDataException("Ungültige SHA256SUMS-Zeile im Runtime-Paket.");
            if (!checksums.TryAdd(match.Groups["name"].Value, match.Groups["hash"].Value))
                throw new InvalidDataException("Doppelter Dateieintrag in Runtime-SHA256SUMS.");
        }
        foreach (var path in new[] { archive, wheel })
        {
            var name = Path.GetFileName(path);
            if (!checksums.TryGetValue(name, out var expected))
                throw new InvalidDataException($"SHA256-Nachweis fehlt für {name}.");
            await using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read,
                1024 * 1024, useAsync: true);
            var actual = Convert.ToHexString(await SHA256.HashDataAsync(stream, cancellationToken));
            if (!string.Equals(actual, expected, StringComparison.OrdinalIgnoreCase))
                throw new InvalidDataException($"SHA256-Prüfung fehlgeschlagen: {name}. Installation nicht gestartet.");
        }
        return new RuntimePackage(archive, wheel);
    }
}
