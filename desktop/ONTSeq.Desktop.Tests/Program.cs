using ONTSeq.Desktop;
using System.Diagnostics;
using System.Net;
using System.Net.Sockets;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;

const string UnicodeOutputFixture = "DNS-Auflösung fehlgeschlagen – Grüße, Straße, Methylierung 🧬\n";
if (args is ["--emit-utf8-fixture"])
{
    var bytes = Encoding.UTF8.GetBytes(UnicodeOutputFixture);
    await Console.OpenStandardOutput().WriteAsync(bytes);
    await Console.OpenStandardError().WriteAsync(bytes);
    return;
}
if (args is ["--emit-wsl-utf16-fixture"])
{
    var bytes = Encoding.Unicode.GetBytes("Fehler: WSL-Dienst nicht verfügbar. Error code: Wsl/Service/Test\r\n");
    await Console.OpenStandardOutput().WriteAsync(bytes);
    await Console.OpenStandardError().WriteAsync(bytes);
    return;
}
if (args is ["--verify-wsl-unicode", var unicodeDistribution])
{
    var command = WslServiceLauncher.BackendInvocation(new DesktopSettings { BackendCommand = "python3" },
        "-c", "import sys; assert sys.stdout.encoding.lower()=='utf-8'; assert sys.stderr.encoding.lower()=='utf-8'; s='DNS-Auflösung fehlgeschlagen – Grüße, Straße, Methylierung 🧬\\n'; sys.stdout.write(s); sys.stderr.write(s)");
    var output = await WslServiceLauncher.RunWslAsync(unicodeDistribution, command, CancellationToken.None);
    AssertEqual("0", output.ExitCode.ToString(), "actual WSL unicode producer exit code");
    AssertEqual(UnicodeOutputFixture, output.StdOut, "actual WSL stdout preserves Unicode exactly");
    AssertEqual(UnicodeOutputFixture, output.StdErr, "actual WSL stderr preserves Unicode exactly");
    var missing = await WslServiceLauncher.RunWslAsync("ONTSeq-Nonexistent-Unicode-Test-Distribution",
        ["true"], CancellationToken.None);
    AssertEqual("False", (missing.ExitCode == 0).ToString(), "actual WSL infrastructure error exit code");
    AssertEqual("False", (missing.StdOut + missing.StdErr).Contains('\0').ToString(),
        "actual WSL infrastructure output is decoded rather than exposing UTF-16 nulls");
    Console.WriteLine("Actual WSL UTF-8 stdout/stderr and infrastructure diagnostic decoding passed.");
    return;
}

// Packaging verification hook: hashes the real artifacts and emits the same installer
// command used by Desktop, but never runs it or changes Desktop settings.
if (args is ["--emit-verified-runtime-install", var archivePath, var newPrefix])
{
    var verifiedPackage = await RuntimePackage.VerifyAsync(archivePath, "0.7.1", CancellationToken.None);
    Console.WriteLine(WslServiceLauncher.RuntimeInstallCoreCommand(verifiedPackage, newPrefix));
    return;
}
if (args is ["--verify-installed-runtime", var distribution, var installedPrefix])
{
    var result = await WslServiceLauncher.RunWslAsync(distribution,
        ["sh", "-lc", WslServiceLauncher.RuntimeVersionCheckCommand(installedPrefix)],
        CancellationToken.None);
    Console.Write(result.StdOut);
    Console.Error.Write(result.StdErr);
    Environment.ExitCode = result.ExitCode;
    return;
}

var root = Path.Combine(Path.GetTempPath(), "ONTSeq.Desktop.Tests", Guid.NewGuid().ToString("N"));
Directory.CreateDirectory(root);

try
{
    var utf8StartInfo = WslServiceLauncher.WslProcessStartInfo("synthetic", ["printf", "unused"]);
    AssertEqual("utf-8", utf8StartInfo.StandardOutputEncoding?.WebName,
        "Linux stdout uses explicit UTF-8 independently of Windows code page");
    AssertEqual("utf-8", utf8StartInfo.StandardErrorEncoding?.WebName,
        "Linux stderr uses explicit UTF-8 independently of Windows code page");
    foreach (var (hook, expected) in new[] {
        ("--emit-utf8-fixture", UnicodeOutputFixture),
        ("--emit-wsl-utf16-fixture", "Fehler: WSL-Dienst nicht verfügbar. Error code: Wsl/Service/Test\r\n")
    })
    {
        var childInfo = WslServiceLauncher.WslProcessStartInfo("unused", []);
        childInfo.FileName = Environment.ProcessPath ?? throw new InvalidOperationException("Missing test host path");
        childInfo.ArgumentList.Clear();
        if (Path.GetFileNameWithoutExtension(childInfo.FileName).Equals("dotnet", StringComparison.OrdinalIgnoreCase))
            childInfo.ArgumentList.Add(System.Reflection.Assembly.GetExecutingAssembly().Location);
        childInfo.ArgumentList.Add(hook);
        using var child = Process.Start(childInfo) ?? throw new InvalidOperationException("Unable to start Unicode test child");
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(15));
        var stdout = WslServiceLauncher.ReadWslOutputAsync(child.StandardOutput.BaseStream, timeout.Token);
        var stderr = WslServiceLauncher.ReadWslOutputAsync(child.StandardError.BaseStream, timeout.Token);
        await child.WaitForExitAsync(timeout.Token);
        AssertEqual("0", child.ExitCode.ToString(), "native output producer exits normally");
        AssertEqual(expected, await stdout, "native redirected stdout preserves exact Unicode bytes");
        AssertEqual(expected, await stderr, "native redirected stderr preserves exact Unicode bytes");
    }
    foreach (var (encoding, prefix, input) in new (Encoding, byte[], string)[] {
        (new UTF8Encoding(false, true), [], new string('a', 4099) + "ö🧬"),
        (new UnicodeEncoding(false, false, true), [], new string('a', 2049) + "🧬ö"),
        (new UnicodeEncoding(true, true, true), [0xfe, 0xff], UnicodeOutputFixture),
        (new UTF8Encoding(true, true), [0xef, 0xbb, 0xbf], UnicodeOutputFixture)
    })
    {
        using var stream = new MemoryStream(prefix.Concat(encoding.GetBytes(input)).ToArray());
        AssertEqual(input, await WslServiceLauncher.ReadWslOutputAsync(stream, CancellationToken.None),
            "byte-order marks and split multibyte characters decode without replacement");
    }
    AssertSequenceEqual(
        ["env", "PYTHONIOENCODING=utf-8", "PYTHONUTF8=1", "ontseq", "--version"],
        WslServiceLauncher.BackendInvocation(new DesktopSettings(), "--version"),
        "known Python backend receives explicit Unicode producer environment as separate arguments");

    var methylationProbe = new MethylationProbeResponse(
        "/approved/sample.bam", "detected", "Synthetic test probe", 10, false);
    AssertEqual("detected", methylationProbe.RequireBam("/approved/sample.bam").Status,
        "methylation detection remains bound to selected BAM");
    AssertThrows<InvalidDataException>(() => methylationProbe.RequireBam("/approved/other.bam"),
        "a different BAM cannot reuse detection");
    AssertThrows<InvalidDataException>(
        () => (methylationProbe with { Status = "not_detected", Complete = false }).RequireBam("/approved/sample.bam"),
        "a partial negative cannot be called not detected");
    var incompleteProbe = methylationProbe with
    {
        Status = "unknown", ReasonCode = "sample_incomplete", CheckedReads = 10000,
        ElapsedSeconds = 2.5, ScanMode = "quick"
    };
    var incompletePresentation = MethylationProbePresentation.From(incompleteProbe);
    AssertEqual("Methylierungsstatus offen", incompletePresentation.Title,
        "bounded negative remains unknown in the German interface");
    AssertEqual("True", incompletePresentation.Explanation.Contains("keine sichere Aussage", StringComparison.Ordinal).ToString(),
        "unknown does not imply methylation is absent");
    AssertEqual("10.000 Reads geprüft · 2,5 Sekunden", incompletePresentation.Progress,
        "actual counts and duration use readable German formatting");
    AssertEqual("True", incompletePresentation.Detail.Contains("sample_incomplete", StringComparison.Ordinal).ToString(),
        "machine reason is retained separately from German explanation");
    AssertEqual("Keine Methylierungsinformationen erkannt", MethylationProbePresentation.From(
        incompleteProbe with { Status = "not_detected", Complete = true, ReasonCode = "complete_no_tags" }).Title,
        "only a complete negative is presented as no methylation information");
    AssertEqual("5mC-Methylierungsinformationen erkannt", MethylationProbePresentation.From(
        methylationProbe with { ReasonCode = "detected_5mc", MethylationAvailable = true }).Title,
        "supported 5mC detection is presented positively");
    foreach (var code in new[] { "sample_incomplete", "timeout", "cancelled", "reader_unavailable", "input_unavailable",
        "file_changed", "record_limit", "invalid_tags", "unsupported_modification", "read_error", "worker_failed", "busy",
        "service_unavailable", "request_timeout" })
    {
        var presentation = MethylationProbePresentation.From(incompleteProbe with { ReasonCode = code });
        AssertEqual("True", presentation.Explanation.Contains("keine sichere Aussage", StringComparison.Ordinal).ToString(),
            $"{code} preserves uncertainty");
        AssertEqual("False", presentation.Explanation.StartsWith("Die Vorprüfung konnte", StringComparison.Ordinal).ToString(),
            $"{code} has a specific explanation rather than a generic fallback");
    }
    AssertEqual("BAM-Lesewerkzeug nicht verfügbar", MethylationProbePresentation.From(
        incompleteProbe with { ReasonCode = "reader_unavailable" }).Title,
        "missing reader is distinguished from absent tags");
    AssertEqual("BAM-Lesefehler", MethylationProbePresentation.From(
        incompleteProbe with { ReasonCode = "read_error" }).Title,
        "read failure has a distinct fault title");
    foreach (var code in new[] { "wsl_unavailable", "backend_unavailable", "input_path_unsupported",
        "output_path_unsupported", "input_drive_unavailable", "input_root_missing", "input_root_unreadable",
        "output_drive_unavailable", "output_unwritable", "resource_root_missing", "resource_root_unreadable",
        "runtime_file_missing", "runtime_file_unreadable", "runtime_tool_unavailable", "reference_lock_unavailable",
        "target_bed_unavailable", "resource_profile_unavailable", "prerequisite_check_failed" })
    {
        var failure = MethylationProbeFailures.FromException("/synthetic/input.bam",
            new WslPrerequisiteException(code, "Synthetic preparation detail for " + code), "quick");
        var presentation = MethylationProbePresentation.From(failure);
        AssertEqual(code, failure.ReasonCode, "typed prerequisite reason survives UI adaptation");
        AssertEqual("unknown", failure.Status, "prerequisite failure never becomes absent methylation");
        AssertEqual("True", failure.PreparationFailed.ToString(), "typed prerequisite failure is a preparation error");
        AssertEqual("True", presentation.Explanation.StartsWith("Die BAM-Prüfung hat noch nicht begonnen.", StringComparison.Ordinal).ToString(),
            $"{code} explicitly precedes BAM scanning");
        AssertEqual("False", presentation.Title.Contains("fehlgeschlagen", StringComparison.Ordinal).ToString(),
            $"{code} is not presented as a failed scan worker");
        AssertEqual("Prüfung nicht gestartet · noch keine BAM-Reads geprüft", presentation.Progress,
            $"{code} does not fabricate a zero-second scan result");
        AssertEqual("True", presentation.Detail.Contains(failure.Reason, StringComparison.Ordinal).ToString(),
            $"{code} preserves the precise launcher diagnostic");
    }
    var disconnectedInput = MethylationProbePresentation.From(MethylationProbeFailures.FromException(
        "/synthetic/input.bam", new WslPrerequisiteException("input_drive_unavailable", "Synthetic ENODEV"), "quick"));
    AssertEqual("True", disconnectedInput.Explanation.Contains("unterbrochen", StringComparison.Ordinal).ToString(),
        "unavailable drive can be a broken existing mount rather than only an absent mount");
    foreach (var (error, expectedCode) in new (Exception, string)[] {
        (new InvalidOperationException("Synthetic bootstrap failure"), "preparation_failed"),
        (new System.Net.Http.HttpRequestException("Synthetic connection failure"), "service_unavailable"),
        (new TimeoutException("Synthetic bootstrap timeout"), "request_timeout"),
        (new TaskCanceledException("Synthetic request timeout"), "request_timeout"),
        (new MethylationServiceException("Synthetic busy", "busy"), "busy")
    })
    {
        var preparation = MethylationProbeFailures.FromException("/synthetic/input.bam", error, "thorough", beforeProbe: true);
        AssertEqual(expectedCode, preparation.ReasonCode, "pre-service exception keeps specific transport/service codes");
        AssertEqual("True", preparation.PreparationFailed.ToString(), "bootstrap errors precede any BAM probe");
        AssertEqual("Prüfung nicht gestartet · noch keine BAM-Reads geprüft", MethylationProbePresentation.From(preparation).Progress,
            "bootstrap exceptions are preparation errors even when their reason is transport-specific");
    }
    var workerFailure = MethylationProbeFailures.FromException("/synthetic/input.bam",
        new InvalidOperationException("Synthetic worker failure"), "thorough");
    AssertEqual("worker_failed", workerFailure.ReasonCode, "actual worker failure retains its original classification");
    AssertEqual("False", workerFailure.PreparationFailed.ToString(), "actual worker failure is not relabelled as preparation");
    AssertEqual("Prüfung fehlgeschlagen", MethylationProbePresentation.From(workerFailure).Title,
        "worker result presentation remains unchanged");
    var preparationJson = JsonSerializer.Serialize(MethylationProbeFailures.FromException("/synthetic/input.bam",
        new InvalidOperationException("Synthetic preparation failure"), "quick", beforeProbe: true), JsonDefaults.Options);
    AssertEqual("False", preparationJson.Contains("preparationFailed", StringComparison.OrdinalIgnoreCase).ToString(),
        "Desktop preparation phase is not added to the service's probe schema");
    var outputPathAsFile = Path.Combine(root, "synthetic-output-is-a-file");
    File.WriteAllText(outputPathAsFile, "Synthetic file must remain unchanged");
    try
    {
        DesktopOutputDirectory.EnsureExists(outputPathAsFile);
        throw new InvalidOperationException("Expected a typed output-directory failure");
    }
    catch (WslPrerequisiteException error)
    {
        AssertEqual("output_unwritable", error.ReasonCode, "actual file-instead-of-directory failure has a precise code");
        AssertEqual(outputPathAsFile, error.WindowsPath, "output failure carries the exact configured Windows path");
        AssertEqual("True", (error.InnerException is IOException).ToString(), "output preparation retains its original filesystem exception");
        var presentation = MethylationProbePresentation.From(MethylationProbeFailures.FromException(
            "/synthetic/input.bam", error, "quick", beforeProbe: true));
        AssertEqual("Ausgabeordner nicht nutzbar", presentation.Title,
            "actual invalid output directory is not blamed on the BAM scanner");
        AssertEqual("True", presentation.Explanation.Contains("vorhandene Datei", StringComparison.Ordinal).ToString(),
            "output failure tells the user to select a usable directory");
    }
    AssertEqual("Synthetic file must remain unchanged", File.ReadAllText(outputPathAsFile),
        "output preparation never overwrites a colliding file");
    var preparedOutput = Path.Combine(root, "synthetic-new-output");
    DesktopOutputDirectory.EnsureExists(preparedOutput);
    AssertEqual("True", Directory.Exists(preparedOutput).ToString(), "valid output directories can still be prepared");
    var unsupportedRelativeOutput = Path.GetRelativePath(Environment.CurrentDirectory,
        Path.Combine(root, "synthetic-relative-output-must-not-be-created"));
    var unsupportedOutputTarget = Path.GetFullPath(unsupportedRelativeOutput);
    AssertEqual("False", Directory.Exists(unsupportedOutputTarget).ToString(), "unsupported output fixture starts absent");
    try
    {
        DesktopOutputDirectory.EnsureExists(unsupportedRelativeOutput);
        throw new InvalidOperationException("Expected unsupported output path rejection before directory creation");
    }
    catch (WslPrerequisiteException error)
    {
        AssertEqual("output_path_unsupported", error.ReasonCode, "relative output path is rejected before filesystem mutation");
        AssertEqual(unsupportedRelativeOutput, error.WindowsPath, "rejected configured path is preserved for diagnosis");
    }
    AssertEqual("False", Directory.Exists(unsupportedOutputTarget).ToString(),
        "unsupported output path creates no directory before rejection");
    var prerequisiteSettings = new DesktopSettings
    {
        OutputDirectoryWindows = @"C:\Synthetic results\Output",
        ResourceRootWsl = "/mnt/c/Synthetic resources",
        RuntimeBinWsl = "/home/synthetic/Runtime with spaces/bin"
    };
    var prerequisiteChecks = WslServiceLauncher.ProfilePrerequisiteChecks(prerequisiteSettings, @"E:\Synthetic BAM folder");
    WslServiceLauncher.RequirePrerequisites(prerequisiteChecks, 0, "ONTSEQ_PREREQUISITES_OK\n");
    for (var index = 0; index < prerequisiteChecks.Count; index++)
    {
        var expected = prerequisiteChecks[index];
        try
        {
            WslServiceLauncher.RequirePrerequisites(prerequisiteChecks, 42, $"ONTSEQ_PREREQUISITE:{index}\n");
            throw new InvalidOperationException("Expected a precisely attributed prerequisite failure");
        }
        catch (WslPrerequisiteException error)
        {
            AssertEqual(expected.ReasonCode, error.ReasonCode, "every prerequisite index preserves its typed dependency");
            AssertEqual(expected.WslPath, error.WslPath, "every prerequisite index preserves its WSL path");
            AssertEqual(expected.WindowsPath, error.WindowsPath, "every prerequisite index preserves its Windows path");
            AssertEqual(expected.FailureMessage, error.Message, "every prerequisite index preserves its concrete action");
            if (error.ReasonCode.StartsWith("input_", StringComparison.Ordinal))
                AssertEqual("False", error.Message.Contains("installier", StringComparison.OrdinalIgnoreCase).ToString(),
                    "input access failure does not recommend reinstalling reference resources");
        }
    }
    foreach (var (exitCode, output) in new[] { (0, ""), (0, "ONTSEQ_PREREQUISITE:0"),
        (42, "ONTSEQ_PREREQUISITE:-1"), (42, "ONTSEQ_PREREQUISITE:99999"),
        (42, "ONTSEQ_PREREQUISITE:0\nforeign message"), (1, "ONTSEQ_PREREQUISITES_OK") })
    {
        try
        {
            WslServiceLauncher.RequirePrerequisites(prerequisiteChecks, exitCode, output);
            throw new InvalidOperationException("Expected rejection of invalid prerequisite protocol");
        }
        catch (WslPrerequisiteException error)
        {
            AssertEqual("prerequisite_check_failed", error.ReasonCode,
                "ambiguous prerequisite output cannot be attributed to an arbitrary dependency");
        }
    }
    AssertEqual("True", prerequisiteChecks.Any(check => check.ReasonCode == "output_unwritable" &&
        check.WslPath == "/mnt/c/Synthetic results/Output" && check.FailureMessage.Contains(check.WslPath, StringComparison.Ordinal)).ToString(),
        "output access diagnostic names the actual translated output directory");
    foreach (var invalidPath in new[] { "relative/input", @"C:relative", "/mnt/e/already-translated", @"\rooted", @"\\server\share\input", "C:\\invalid\ninput" })
        AssertThrows<InvalidOperationException>(() => PathBridge.WindowsToWsl(invalidPath),
            "relative, Linux, UNC and control-character paths are rejected without inspection");
    AssertEqual("/mnt/e/Synthetic BAM folder", PathBridge.WindowsToWsl(@"E:\Synthetic BAM folder"),
        "supported absolute drive paths preserve spaces without filesystem lookup");
    var probeGuard = new MethylationProbeGuard();
    var oldScope = probeGuard.Capture("/synthetic/same.bam", "profile-one", "service-one");
    AssertEqual("True", probeGuard.IsCurrent(oldScope, oldScope.BamPath, oldScope.ProfileId, oldScope.ServiceInstanceId).ToString(),
        "current probe scope accepts its own response");
    probeGuard.Invalidate();
    AssertEqual("False", probeGuard.IsCurrent(oldScope, oldScope.BamPath, oldScope.ProfileId, oldScope.ServiceInstanceId).ToString(),
        "reselecting the same BAM or cancelling rejects a late response");
    var currentScope = probeGuard.Capture(oldScope.BamPath, oldScope.ProfileId, oldScope.ServiceInstanceId);
    AssertEqual("False", probeGuard.IsCurrent(currentScope, currentScope.BamPath, "different-profile", currentScope.ServiceInstanceId).ToString(),
        "profile changes invalidate response scope");
    AssertEqual("False", probeGuard.IsCurrent(currentScope, currentScope.BamPath, currentScope.ProfileId, "replacement-service").ToString(),
        "replacement service cannot supply an old selection response");
    var scan = new MethylationScanSnapshot(new string('a', 32), "/approved/sample.bam", "running", 1024, 1.25, null);
    AssertEqual("running", scan.RequireScan(scan.BamPath, scan.ScanId).State, "typed running progress accepted");
    AssertThrows<InvalidDataException>(() => scan.RequireScan("/approved/other.bam", scan.ScanId),
        "scan progress from another BAM is rejected");
    AssertThrows<InvalidDataException>(() => scan.RequireScan(scan.BamPath, new string('b', 32)),
        "scan progress from another job is rejected");
    AssertThrows<InvalidDataException>(() => (scan with { State = "completed" }).RequireScan(scan.BamPath),
        "completed scan must contain a typed result");
    AssertThrows<InvalidDataException>(() => (scan with { ElapsedSeconds = double.NaN }).RequireScan(scan.BamPath),
        "invalid progress duration is rejected");
    AssertThrows<InvalidDataException>(() => MethylationScanSnapshot.RequireId("../other"),
        "scan IDs cannot escape the fixed API route");
    AssertEqual("cancelled", (scan with { State = "cancelled", Result = incompleteProbe with { ReasonCode = "cancelled" } })
        .RequireScan(scan.BamPath).State, "cancelled worker remains unknown rather than becoming a negative result");

    await VerifyMethylationHttpContractAsync();
    VerifyMethylationLayout();
    var methylationRequest = new RunStartRequest("/approved/sample.bam", "sample-001", null,
        "AML_LCWGS_GRCh38", "GRCh38", "lcwgs");
    using (var defaultRequest = JsonDocument.Parse(JsonSerializer.Serialize(methylationRequest, JsonDefaults.Options)))
        AssertEqual("False", defaultRequest.RootElement.GetProperty("include_methylation").GetBoolean().ToString(),
            "methylation defaults to explicit false");
    using (var optedInRequest = JsonDocument.Parse(JsonSerializer.Serialize(
        methylationRequest with { IncludeMethylation = true }, JsonDefaults.Options)))
        AssertEqual("True", optedInRequest.RootElement.GetProperty("include_methylation").GetBoolean().ToString(),
            "explicit methylation decision reaches the backend contract");

    var bam = Path.Combine(root, "sample.bam");
    var preferred = bam + ".bai";
    var alternative = Path.ChangeExtension(bam, ".bai");
    File.WriteAllText(bam, "BAM");

    AssertEqual(null, BamIndexLocator.Find(bam), "missing index");

    File.WriteAllText(alternative, "BAI");
    AssertEqual(alternative, BamIndexLocator.Find(bam), "sample.bai");

    File.WriteAllText(preferred, "BAI");
    AssertEqual(preferred, BamIndexLocator.Find(bam), "sample.bam.bai precedence");

    File.Delete(preferred);
    AssertEqual(alternative, BamIndexLocator.Find(bam), "fallback after preferred removal");

    File.Delete(alternative);
    File.WriteAllText(Path.Combine(root, "other.bai"), "BAI");
    AssertEqual(null, BamIndexLocator.Find(bam), "unrelated index");

    AssertEqual(
        "/mnt/c/Lab/sample.bam",
        PathBridge.WindowsToWsl(@"C:\Lab\sample.bam"),
        "drive path translation");
    AssertEqual(
        "/mnt/p/Lab FG06/NANOPORE/sample.bam",
        PathBridge.WindowsToWsl(@"P:\Lab FG06\NANOPORE\sample.bam"),
        "mapped drive translation");
    AssertThrows<InvalidOperationException>(
        () => PathBridge.WindowsToWsl(@"\\server\share\sample.bam"),
        "UNC refusal");

    var profileDefaults = new DesktopSettings();
    AssertEqual(
        "~/.local/share/ontseq/resources-v" + DesktopVersion.Value.Replace("-engineering", ""),
        profileDefaults.ResourceRootWsl,
        "user-writable default WSL resource root follows the software release");
    AssertEqual(profileDefaults.ResourceRootWsl, DesktopSettings.NormalizeResourceRootWsl(null),
        "missing resource root uses the current release");
    AssertEqual(profileDefaults.ResourceRootWsl, DesktopSettings.NormalizeResourceRootWsl("  "),
        "blank resource root uses the current release");
    var savedResourceSettings = JsonSerializer.Deserialize<DesktopSettings>(
        """{"resourceRootWsl":"/home/synthetic/.local/share/ontseq/resources-v0.6.2-hg19-as"}""",
        JsonDefaults.Options)!;
    savedResourceSettings.ApplyProfileDefaults();
    AssertEqual("/home/synthetic/.local/share/ontseq/resources-v0.6.2-hg19-as",
        savedResourceSettings.ResourceRootWsl,
        "an explicit existing resource installation is not silently relocated");
    AssertEqual("AML_LCWGS_GRCh38", profileDefaults.DefaultProfile, "default analysis profile");
    AssertEqual("8", DesktopProfiles.Supported.Count.ToString(), "exact supported profile count");
    AssertEqual(
        "4",
        DesktopProfiles.Supported.Count(profile => profile.GenomeBuild == "GRCh38").ToString(),
        "four existing GRCh38 profiles remain unchanged");
    AssertSequenceEqual(
        [
            DesktopProfiles.DefaultProfileId,
            DesktopProfiles.AdaptiveSamplingProfileId,
            DesktopProfiles.Canonical25LcwgsProfileId,
            DesktopProfiles.Canonical25AdaptiveSamplingProfileId
        ],
        DesktopProfiles.Supported
            .Where(profile => profile.GenomeBuild == "GRCh38")
            .Select(profile => profile.ProfileId)
            .ToArray(),
        "GRCh38 profile identities and ordering remain unchanged");
    AssertEqual(
        "4",
        DesktopProfiles.Supported.Count(profile => profile.GenomeBuild == "GRCh37").ToString(),
        "GRCh37 exposes lcWGS and Adaptive Sampling for both dictionary contracts");
    var native37Profile = DesktopProfiles.Require(DesktopProfiles.Grch37LcwgsProfileId);
    AssertEqual("GRCh37", native37Profile.GenomeBuild, "native37 build identity");
    AssertEqual("exact_full", native37Profile.DictionaryContract, "native37 exact full assembly contract");
    AssertEqual("lcwgs", native37Profile.Assay, "native37 lcWGS assay");
    AssertEqual("False", native37Profile.AdaptiveSampling.ToString(), "native37 lcWGS remains panel-free");
    AssertEqual("True", native37Profile.DictionaryLabel.Contains("Patches/Haplotypen").ToString(),
        "native37 label does not incorrectly claim primary-only assembly");
    var hg19Profile = DesktopProfiles.Require(
        DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId);
    AssertEqual("GRCh37", hg19Profile.GenomeBuild, "hg19 Canonical-25 build identity");
    AssertEqual("grch37_ucsc_hg19_canonical_25", hg19Profile.DictionaryContract,
        "hg19 profile exposes its strict dictionary contract");
    AssertEqual("True", hg19Profile.DictionaryLabel.Contains("chrM=16571").ToString(),
        "hg19 profile makes the mitochondrial length visible");
    AssertEqual(
        "AML_AS_111_GRCh38",
        DesktopProfiles.Require("AML_AS_111_GRCh38").ProfileId,
        "adaptive sampling profile identity");
    AssertEqual(
        "Canonical-25 (chr1–22, chrX, chrY, chrM)",
        DesktopProfiles.Require("AML_LCWGS_GRCh38_CANONICAL25").DictionaryLabel,
        "Canonical-25 profile exposes its exact BAM dictionary contract");
    var native37AdaptiveProfile = DesktopProfiles.Require(
        DesktopProfiles.Grch37AdaptiveSamplingProfileId);
    AssertEqual("GRCh37", native37AdaptiveProfile.GenomeBuild,
        "native37 Adaptive Sampling build identity");
    AssertEqual("adaptive_sampling", native37AdaptiveProfile.Assay,
        "native37 Adaptive Sampling assay binding");
    AssertEqual("exact_full", native37AdaptiveProfile.DictionaryContract,
        "native37 Adaptive Sampling full publisher dictionary contract");
    AssertEqual("True", native37AdaptiveProfile.AdaptiveSampling.ToString(),
        "native37 Adaptive Sampling capability flag");
    AssertEqual("True", native37AdaptiveProfile.DisplayName.Contains(
        "110/111 kartiert", StringComparison.Ordinal).ToString(),
        "native37 Adaptive Sampling label discloses mapped target count");
    var hg19AdaptiveProfile = DesktopProfiles.Require(
        DesktopProfiles.Grch37UcscHg19Canonical25AdaptiveSamplingProfileId);
    AssertEqual("GRCh37", hg19AdaptiveProfile.GenomeBuild,
        "hg19 Adaptive Sampling build identity");
    AssertEqual("adaptive_sampling", hg19AdaptiveProfile.Assay,
        "hg19 Adaptive Sampling assay binding");
    AssertEqual("grch37_ucsc_hg19_canonical_25", hg19AdaptiveProfile.DictionaryContract,
        "hg19 Adaptive Sampling strict dictionary contract");
    AssertEqual("True", hg19AdaptiveProfile.DictionaryLabel.Contains(
        "chrM=16571", StringComparison.Ordinal).ToString(),
        "hg19 Adaptive Sampling makes the mitochondrial dictionary identity visible");
    AssertEqual("True", hg19AdaptiveProfile.DisplayName.Contains(
        "110/111 kartiert", StringComparison.Ordinal).ToString(),
        "hg19 Adaptive Sampling label discloses mapped target count");
    AssertSequenceEqual(
        [
            "GRCh38_GENCODE50_MANE1.5_v1",
            "HEMATOLOGY_v3",
            "AML_AS_111_GRCh38_v1"
        ],
        WslServiceLauncher.ManagedGrch38ResourceBundleIds,
        "Desktop repair owns the complete GRCh38 profile resource family");
    AssertSequenceEqual(
        ["GRCh37_GENCODE19_HG19_v2", "HEMATOLOGY_GRCh37_v1", "AML_AS_111_GRCh37_v1"],
        WslServiceLauncher.ManagedResourceBundleIds("GRCh37"),
        "native37 family includes its build-specific Adaptive Sampling panel");
    AssertEqual("False", WslServiceLauncher.ManagedResourceBundleIds("GRCh37").Contains(
        WslServiceLauncher.AmlAdaptivePanelBundleId, StringComparer.Ordinal).ToString(),
        "native37 resources never reuse the GRCh38 Adaptive Sampling panel");
    AssertThrows<ArgumentOutOfRangeException>(
        () => WslServiceLauncher.ManagedResourceBundleIds("hg19"), "no implicit build alias fallback");

    var installedRuntimeSettings = new DesktopSettings
    {
        RuntimeBinWsl = "/opt/ontseq/bin"
    };
    const string runtimeShare = "/opt/ontseq/share/ontseq/";
    AssertSequenceEqual(
        [
            runtimeShare + "configs/qc/defaults.yaml",
            runtimeShare + "configs/qc/adaptive_target_coverage.technical.yaml",
            runtimeShare + "configs/components/default.yaml",
            runtimeShare + "configs/sv/sniffles2.conservative.technical.yaml",
            runtimeShare + "configs/sv/cutesv.conservative.technical.yaml",
            runtimeShare + "configs/sv/sniffles2_cutesv.consensus.technical.yaml",
            runtimeShare + "configs/sv/evidence-priority.technical.yaml",
            runtimeShare + "configs/cnv/qdnaseq_ace.technical.yaml",
            runtimeShare + "scripts/run_qdnaseq_ace.R"
        ],
        WslServiceLauncher.RequiredRuntimeFiles(installedRuntimeSettings),
        "complete packaged runtime file contract");
    AssertSequenceEqual(
        [
            "/opt/ontseq/bin/ontseq",
            "/opt/ontseq/bin/Rscript",
            "/opt/ontseq/bin/samtools",
            "/opt/ontseq/bin/cramino",
            "/opt/ontseq/bin/sniffles",
            "/opt/ontseq/bin/cuteSV",
            "/opt/ontseq/bin/mosdepth"
        ],
        WslServiceLauncher.RequiredRuntimeTools(installedRuntimeSettings),
        "complete packaged runtime tool contract");
    AssertSequenceEqual(
        [
            "--qc-policy", runtimeShare + "configs/qc/defaults.yaml",
            "--sniffles-policy", runtimeShare + "configs/sv/sniffles2.conservative.technical.yaml",
            "--cutesv-policy", runtimeShare + "configs/sv/cutesv.conservative.technical.yaml",
            "--sv-consensus-policy", runtimeShare + "configs/sv/sniffles2_cutesv.consensus.technical.yaml",
            "--sv-evidence-policy", runtimeShare + "configs/sv/evidence-priority.technical.yaml",
            "--target-coverage-policy", runtimeShare + "configs/qc/adaptive_target_coverage.technical.yaml",
            "--components", runtimeShare + "configs/components/default.yaml",
            "--cnv-policy", runtimeShare + "configs/cnv/qdnaseq_ace.technical.yaml",
            "--qdnaseq-rscript", "/opt/ontseq/bin/Rscript",
            "--qdnaseq-script", runtimeShare + "scripts/run_qdnaseq_ace.R"
        ],
        WslServiceLauncher.BundledPolicyArguments(
            installedRuntimeSettings,
            includeCnv: true,
            includeCore034: true),
        "profile service receives absolute packaged policy paths");
    AssertSequenceEqual(
        [
            "--qc-policy", runtimeShare + "configs/qc/defaults.yaml",
            "--sniffles-policy", runtimeShare + "configs/sv/sniffles2.conservative.technical.yaml",
            "--cnv-policy", runtimeShare + "configs/cnv/qdnaseq_ace.technical.yaml",
            "--qdnaseq-rscript", "/opt/ontseq/bin/Rscript",
            "--qdnaseq-script", runtimeShare + "scripts/run_qdnaseq_ace.R"
        ],
        WslServiceLauncher.BundledPolicyArguments(
            installedRuntimeSettings,
            includeCnv: true,
            includeCore034: false),
        "system smoke retains its smaller accepted argument contract");
    AssertSequenceEqual(
        [],
        WslServiceLauncher.BundledPolicyArguments(
            new DesktopSettings(),
            includeCnv: true,
            includeCore034: true),
        "external development backend retains its own policy defaults");
    AssertEqual(
        runtimeShare + "configs/qc/defaults.yaml",
        WslServiceLauncher.RequiredRuntimeFiles(
            new DesktopSettings { RuntimeBinWsl = "/opt/ontseq/bin/" })[0],
        "runtime bin trailing slash does not move the packaged asset root");
    AssertEqual(
        "/srv/ontseq",
        DesktopSettings.NormalizeResourceRootWsl(" /srv/ontseq/ "),
        "resource root normalization");
    AssertEqual(
        "~/.local/share/ontseq/resources",
        DesktopSettings.NormalizeResourceRootWsl(" ~/.local/share/ontseq/resources/ "),
        "home-relative resource root normalization");
    AssertThrows<InvalidDataException>(
        () => DesktopSettings.NormalizeResourceRootWsl("relative/resources"),
        "relative resource root refusal");
    AssertThrows<InvalidDataException>(
        () => DesktopSettings.NormalizeResourceRootWsl("/opt/../mixed-build-root"),
        "relative WSL segment refusal");

    const string separateModkitPath = "/opt/local tools/modkit 0.4.1/bin/modkit";
    AssertEqual(null, DesktopSettings.ValidateModkitExecutableWsl(null),
        "modkit null preserves PATH lookup");
    AssertEqual(separateModkitPath, DesktopSettings.ValidateModkitExecutableWsl(separateModkitPath),
        "separate modkit executable preserves spaces as an opaque path");
    foreach (var unsafeModkitPath in new[] {
        "", "/", "/opt/tools/", "modkit", "~/bin/modkit", "../modkit",
        "/opt/../bin/modkit", "/opt/./bin/modkit", "//server/modkit",
        @"C:\tools\modkit.exe", "/opt/tools\n/modkit", "/opt/tools\t/modkit",
        "/opt/$(touch bad)/modkit", "/opt/modkit;echo bad", "/opt/modkit|cat",
        "/opt/`id`/modkit", "/opt/modkit && id", "/opt/modkit>bad"
    })
    {
        AssertThrows<InvalidDataException>(
            () => DesktopSettings.ValidateModkitExecutableWsl(unsafeModkitPath),
            "unsafe modkit executable path is refused");
        AssertThrows<InvalidDataException>(
            () => WslServiceLauncher.ProfileServiceArguments(
                new DesktopSettings { ModkitExecutableWsl = unsafeModkitPath }, root, 8765,
                "11111111111111111111111111111111"),
            "direct profile service argument construction cannot bypass modkit validation");
    }
    foreach (var runtimeBin in new string?[] { null, "/opt/ontseq/bin" })
    {
        var separateToolSettings = new DesktopSettings
        {
            RuntimeBinWsl = runtimeBin,
            ModkitExecutableWsl = separateModkitPath,
            OutputDirectoryWindows = Path.Combine(root, "separate-tool-results")
        };
        var profileArguments = WslServiceLauncher.ProfileServiceArguments(
            separateToolSettings, root, 8765, "11111111111111111111111111111111").ToArray();
        AssertEqual("1", profileArguments.Count(item => item == "--modkit").ToString(),
            "profile service receives one explicit modkit argument with or without a bundled runtime");
        AssertEqual(separateModkitPath, profileArguments[Array.IndexOf(profileArguments, "--modkit") + 1],
            "modkit path remains one complete profile argument");
        var legacyArguments = WslServiceLauncher.LegacyServiceArguments(
            separateToolSettings, root, "/synthetic/reference.lock.json").ToArray();
        AssertEqual(separateModkitPath, legacyArguments[Array.IndexOf(legacyArguments, "--modkit") + 1],
            "legacy service start uses the same explicit modkit path");
        AssertEqual("False", WslServiceLauncher.BundledPolicyArguments(
            separateToolSettings, true, false).Contains("--modkit").ToString(),
            "system-smoke policies never acquire an unsupported modkit flag");
    }
    AssertEqual("False", WslServiceLauncher.ProfileServiceArguments(
        new DesktopSettings(), root, 8765, "11111111111111111111111111111111")
        .Contains("--modkit").ToString(), "default service keeps backend modkit PATH behavior");

    AssertSequenceEqual(
        ["references", "status", "--resource-root", "/opt/ontseq"],
        WslServiceLauncher.ResourceManagementArguments("status", "/opt/ontseq/"),
        "bundle status command bridge");
    AssertSequenceEqual(
        [
            "references", "status", "--resource-root",
            profileDefaults.ResourceRootWsl
        ],
        WslServiceLauncher.ResourceManagementArguments(
            "status", DesktopSettings.DefaultResourceRootWsl),
        "home-relative Desktop bundle status command bridge");
    AssertSequenceEqual(
        ["references", "validate", "--resource-root", "/opt/ontseq"],
        WslServiceLauncher.ResourceManagementArguments("validate", "/opt/ontseq"),
        "bundle validation command bridge");
    AssertSequenceEqual(
        [
            "references", "install", "GRCh38_GENCODE50_MANE1.5_v1",
            "--resource-root", "/opt/ontseq"
        ],
        WslServiceLauncher.ResourceManagementArguments(
            "install", "/opt/ontseq", "GRCh38_GENCODE50_MANE1.5_v1"),
        "bundle install command bridge");
    AssertSequenceEqual(
        [
            "references", "repair", "GRCh38_GENCODE50_MANE1.5_v1",
            "--resource-root", "/opt/ontseq"
        ],
        WslServiceLauncher.ResourceManagementArguments(
            "repair", "/opt/ontseq", "GRCh38_GENCODE50_MANE1.5_v1"),
        "bundle repair command bridge");
    foreach (var action in new[] { "install", "repair" })
        AssertSequenceEqual(
            ["references", action, "GRCh37_GENCODE19_HG19_v2", "--resource-root", "/opt/ontseq"],
            WslServiceLauncher.ResourceManagementArguments(action, "/opt/ontseq",
                WslServiceLauncher.ManagedResourceBundleIds("GRCh37")[0]),
            "native37 " + action + " command bridge");

    const string readyResourceStatus = """
        {
          "references": [
            {"bundle_id": "GRCh38_GENCODE50_MANE1.5_v1", "valid": true}
          ],
          "profiles": [
            "AML_LCWGS_GRCh38",
            "AML_AS_111_GRCh38",
            "AML_LCWGS_GRCh38_CANONICAL25",
            "AML_AS_111_GRCh38_CANONICAL25"
          ],
          "diagnostics": []
        }
        """;
    var readyStatus = WslServiceLauncher.InterpretResourceStatus(0, readyResourceStatus, "");
    AssertEqual("True", readyStatus.Ok.ToString(), "ready GRCh38 resource status");
    AssertEqual(
        "True",
        WslServiceLauncher.ManagedGrch38ResourceBundleIds.All(
            bundle => readyStatus.Detail.Contains(bundle, StringComparison.Ordinal)).ToString(),
        "ready status names every repair-managed bundle");
    const string incompleteResourceStatus = """
        {
          "references": [
            {"bundle_id": "GRCh38_GENCODE50_MANE1.5_v1", "valid": true}
          ],
          "profiles": [
            "AML_LCWGS_GRCh38",
            "AML_AS_111_GRCh38",
            "AML_LCWGS_GRCh38_CANONICAL25"
          ],
          "diagnostics": []
        }
        """;
    var incompleteStatus = WslServiceLauncher.InterpretResourceStatus(
        0, incompleteResourceStatus, "");
    AssertEqual("False", incompleteStatus.Ok.ToString(), "missing Canonical-25 AS profile status refusal");
    AssertEqual(
        "True",
        incompleteStatus.Detail.Contains(
            "AML_AS_111_GRCh38_CANONICAL25", StringComparison.Ordinal).ToString(),
        "missing Canonical-25 AS profile is named");

    const string ready37Status = """
        {"references":[{"bundle_id":"GRCh37_GENCODE19_HG19_v2","valid":true}],
         "profiles":["AML_LCWGS_GRCh37","AML_AS_111_GRCh37",
          "AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25",
          "AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25"],"diagnostics":[]}
        """;
    var ready37StatusResult = WslServiceLauncher.InterpretResourceStatus(
        0, ready37Status, "", "GRCh37");
    AssertEqual("True", ready37StatusResult.Ok.ToString(),
        "selected native37 family ready without any GRCh38 resources");
    AssertEqual("True", WslServiceLauncher.ManagedGrch37ResourceBundleIds.All(
        bundle => ready37StatusResult.Detail.Contains(bundle, StringComparison.Ordinal)).ToString(),
        "ready GRCh37 status names reference, knowledge and its own Adaptive Sampling panel");
    AssertEqual("True", WslServiceLauncher.InterpretResourceStatus(0, readyResourceStatus, "", "GRCh38").Ok.ToString(),
        "selected GRCh38 family ready without native37 resources");
    AssertEqual("False", WslServiceLauncher.InterpretResourceStatus(0, readyResourceStatus, "", "GRCh37").Ok.ToString(),
        "GRCh38 resources cannot satisfy native37 readiness");
    AssertEqual("False", WslServiceLauncher.InterpretResourceStatus(0, ready37Status, "", "GRCh38").Ok.ToString(),
        "native37 resources cannot satisfy GRCh38 readiness");
    AssertEqual("False", WslServiceLauncher.InterpretResourceStatus(0,
        ready37Status.Replace("\"valid\":true", "\"valid\":false"), "", "GRCh37").Ok.ToString(),
        "invalid native37 reference rejected");
    AssertEqual("False", WslServiceLauncher.InterpretResourceStatus(0,
        ready37Status.Replace("\"AML_LCWGS_GRCh37\"", "\"AML_LCWGS_GRCh38\""), "", "GRCh37").Ok.ToString(),
        "missing native37 profile rejected");
    AssertEqual("False", WslServiceLauncher.InterpretResourceStatus(0,
        ready37Status.Replace("\"valid\":true", "\"valid\":true,\"genome_build\":\"GRCh38\""), "", "GRCh37").Ok.ToString(),
        "explicit cross-build reference metadata rejected");

    var only37Families = WslServiceLauncher.InterpretResourceFamilyStates(0, ready37Status, "");
    var ready37Family = only37Families["GRCh37"];
    var absent38Family = only37Families["GRCh38"];
    AssertEqual(ResourceFamilyAvailability.Ready.ToString(), ready37Family.Availability.ToString(),
        "only-GRCh37 status classifies GRCh37 as ready");
    AssertEqual("True", ready37Family.CanAnalyze.ToString(), "ready family enables analysis");
    AssertEqual("False", ready37Family.CanInstall.ToString(), "ready family disables install");
    AssertEqual("True", ready37Family.CanRepair.ToString(), "ready family permits explicit full repair");
    AssertEqual(ResourceFamilyAvailability.NotInstalled.ToString(), absent38Family.Availability.ToString(),
        "only-GRCh37 status classifies GRCh38 as not installed");
    AssertEqual("False", absent38Family.CanAnalyze.ToString(), "missing family blocks analysis");
    AssertEqual("True", absent38Family.CanInstall.ToString(), "missing family enables install");
    AssertEqual("False", absent38Family.CanRepair.ToString(), "missing family disables repair");

    var only38Families = WslServiceLauncher.InterpretResourceFamilyStates(0, readyResourceStatus, "");
    AssertEqual(ResourceFamilyAvailability.Ready.ToString(),
        only38Families["GRCh38"].Availability.ToString(), "only-GRCh38 status classifies GRCh38 as ready");
    AssertEqual(ResourceFamilyAvailability.NotInstalled.ToString(),
        only38Families["GRCh37"].Availability.ToString(), "only-GRCh38 status leaves GRCh37 uninstalled");
    const string nothingInstalledStatus =
        "{\"references\":[],\"profiles\":[],\"profile_status\":[],\"diagnostics\":[]}";
    var nothingInstalled = WslServiceLauncher.InterpretResourceFamilyStates(
        0, nothingInstalledStatus, "");
    AssertEqual("True", DesktopResourcePolicy.GenomeBuilds.All(build =>
            nothingInstalled[build].Availability == ResourceFamilyAvailability.NotInstalled).ToString(),
        "empty but valid status classifies both families as not installed");

    const string bothFamiliesReadyStatus = """
        {"references":[
           {"bundle_id":"GRCh37_GENCODE19_HG19_v2","valid":true,"genome_build":"GRCh37"},
           {"bundle_id":"GRCh38_GENCODE50_MANE1.5_v1","valid":true,"genome_build":"GRCh38"}],
         "profiles":[
           "AML_LCWGS_GRCh37","AML_AS_111_GRCh37",
           "AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25",
           "AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25",
           "AML_LCWGS_GRCh38","AML_AS_111_GRCh38",
           "AML_LCWGS_GRCh38_CANONICAL25","AML_AS_111_GRCh38_CANONICAL25"],
         "profile_status":[],"diagnostics":[]}
        """;
    var bothFamilies = WslServiceLauncher.InterpretResourceFamilyStates(
        0, bothFamiliesReadyStatus, "");
    AssertEqual(ResourceFamilyAvailability.Ready.ToString(),
        bothFamilies["GRCh37"].Availability.ToString(), "both-ready status keeps GRCh37 ready");
    AssertEqual(ResourceFamilyAvailability.Ready.ToString(),
        bothFamilies["GRCh38"].Availability.ToString(), "both-ready status keeps GRCh38 ready");

    const string incomplete37Status = """
        {"references":[{"bundle_id":"GRCh37_GENCODE19_HG19_v2","valid":true}],
         "profiles":["AML_LCWGS_GRCh37"],"profile_status":[],"diagnostics":[]}
        """;
    var incomplete37 = WslServiceLauncher.InterpretResourceFamilyStates(
        0, incomplete37Status, "")["GRCh37"];
    AssertEqual(ResourceFamilyAvailability.Incomplete.ToString(), incomplete37.Availability.ToString(),
        "missing family profile is incomplete rather than not installed");
    AssertEqual("True", incomplete37.Detail.Contains(
        DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId, StringComparison.Ordinal).ToString(),
        "incomplete status names missing hg19 profile");
    AssertEqual("True", incomplete37.Detail.Contains(
        DesktopProfiles.Grch37AdaptiveSamplingProfileId, StringComparison.Ordinal).ToString(),
        "incomplete status names missing native37 Adaptive Sampling profile");
    AssertEqual("False", incomplete37.CanInstall.ToString(), "incomplete family disables install");
    AssertEqual("True", incomplete37.CanRepair.ToString(), "incomplete family enables repair");

    const string invalid37ReferenceStatus = """
        {"references":[{"bundle_id":"GRCh37_GENCODE19_HG19_v2","valid":false}],
         "profiles":["AML_LCWGS_GRCh37","AML_AS_111_GRCh37",
          "AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25",
          "AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25"],
         "profile_status":[],"diagnostics":[]}
        """;
    var invalid37Reference = WslServiceLauncher.InterpretResourceFamilyStates(
        0, invalid37ReferenceStatus, "")["GRCh37"];
    AssertEqual(ResourceFamilyAvailability.Incomplete.ToString(),
        invalid37Reference.Availability.ToString(), "invalid installed reference requires repair");
    AssertEqual("True", invalid37Reference.Detail.Contains(
        WslServiceLauncher.Grch37ReferenceBundleId, StringComparison.Ordinal).ToString(),
        "invalid installed reference is named");
    var wrongBuild37Reference = WslServiceLauncher.InterpretResourceFamilyStates(
        0,
        ready37Status.Replace(
            "\"valid\":true",
            "\"valid\":true,\"genome_build\":\"GRCh38\"",
            StringComparison.Ordinal),
        "")["GRCh37"];
    AssertEqual(ResourceFamilyAvailability.Incomplete.ToString(),
        wrongBuild37Reference.Availability.ToString(),
        "cross-build metadata is incomplete rather than absent");

    const string orphaned37ProfileStatus = """
        {"references":[],"profiles":[],
         "profile_status":[{"profile_id":"AML_LCWGS_GRCh37","valid":false}],"diagnostics":[]}
        """;
    AssertEqual(ResourceFamilyAvailability.Incomplete.ToString(),
        WslServiceLauncher.InterpretResourceFamilyStates(
            0, orphaned37ProfileStatus, "")["GRCh37"].Availability.ToString(),
        "observed damaged family profile is not misclassified as uninstalled");

    var unavailableFamilies = WslServiceLauncher.InterpretResourceFamilyStates(
        1, "", "backend unavailable");
    foreach (var build in DesktopResourcePolicy.GenomeBuilds)
    {
        var unavailable = unavailableFamilies[build];
        AssertEqual(ResourceFamilyAvailability.Unavailable.ToString(), unavailable.Availability.ToString(),
            build + " command failure status");
        AssertEqual("False", unavailable.CanAnalyze.ToString(), build + " unavailable analysis refusal");
        AssertEqual("False", unavailable.CanInstall.ToString(), build + " unavailable install refusal");
        AssertEqual("False", unavailable.CanRepair.ToString(), build + " unavailable repair refusal");
        AssertEqual("True", unavailable.Detail.Contains("backend unavailable", StringComparison.Ordinal).ToString(),
            build + " command failure detail");
    }
    AssertEqual(ResourceFamilyAvailability.Unavailable.ToString(),
        WslServiceLauncher.InterpretResourceFamilyStates(0, "not-json", "")["GRCh37"].Availability.ToString(),
        "malformed resource status fails closed");
    AssertEqual(ResourceFamilyAvailability.Unavailable.ToString(),
        WslServiceLauncher.InterpretResourceFamilyStates(0, "{}", "")["GRCh38"].Availability.ToString(),
        "missing resource status arrays fail closed");

    var only37ProfileOptions = DesktopResourcePolicy.BuildProfileOptions(only37Families);
    AssertEqual("8", only37ProfileOptions.Count.ToString(), "all profile capabilities remain visible");
    AssertEqual("4", only37ProfileOptions.Count(option => option.IsEnabled).ToString(),
        "only locally ready GRCh37 profiles are enabled");
    AssertEqual("True", only37ProfileOptions
        .Where(option => option.Profile.GenomeBuild == "GRCh37" && option.Profile.AdaptiveSampling)
        .All(option => option.IsEnabled).ToString(),
        "both locally ready GRCh37 Adaptive Sampling contracts are selectable");
    AssertEqual("True", only37ProfileOptions
        .Where(option => option.Profile.GenomeBuild == "GRCh38")
        .All(option => !option.IsEnabled && option.ToString().Contains("nicht installiert", StringComparison.Ordinal))
        .ToString(), "uninstalled GRCh38 choices remain visible and labelled");
    var retainedHg19 = DesktopResourcePolicy.ResolveInitialProfile(
        DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId, only37Families);
    AssertEqual(DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId,
        retainedHg19.Profile?.ProfileId, "ready configured hg19 profile is retained");
    AssertEqual("False", retainedHg19.IsFallback.ToString(), "ready configured profile needs no fallback");
    var fallbackFrom38 = DesktopResourcePolicy.ResolveInitialProfile(
        DesktopProfiles.DefaultProfileId, only37Families);
    AssertEqual(DesktopProfiles.Grch37LcwgsProfileId, fallbackFrom38.Profile?.ProfileId,
        "unavailable configured profile falls back deterministically to first ready profile");
    AssertEqual("True", fallbackFrom38.IsFallback.ToString(), "unavailable configured profile reports fallback");
    AssertEqual("True", fallbackFrom38.Notice?.Contains(
        DesktopProfiles.DefaultProfileId, StringComparison.Ordinal).ToString(),
        "fallback notice names unavailable configured profile");
    AssertEqual("True", DesktopResourcePolicy.CanPersistProfile(
        DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId, only37Families).ToString(),
        "ready profile may be persisted");
    AssertEqual("False", DesktopResourcePolicy.CanPersistProfile(
        DesktopProfiles.DefaultProfileId, only37Families).ToString(),
        "unavailable profile may not be persisted");
    var noReadyDecision = DesktopResourcePolicy.ResolveInitialProfile(
        DesktopProfiles.DefaultProfileId,
        DesktopResourcePolicy.UnavailableFamilies("status unavailable"));
    AssertEqual(null, noReadyDecision.Profile?.ProfileId, "no ready family leaves analysis selection blocked");
    var retainedReady38 = DesktopResourcePolicy.ResolveInitialProfile(
        DesktopProfiles.Canonical25LcwgsProfileId, bothFamilies);
    AssertEqual(DesktopProfiles.Canonical25LcwgsProfileId, retainedReady38.Profile?.ProfileId,
        "both-ready state retains exact configured profile");
    AssertEqual("False", retainedReady38.IsFallback.ToString(), "both-ready exact selection needs no fallback");

    AssertEqual("http://127.0.0.1:8765/workspace", OntSeqServiceClient.WorkspaceUri(8765).AbsoluteUri,
        "live workspace uses loopback service route without token in URL");
    AssertEqual(
        "http://127.0.0.1:18767/workspace?profile=AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25",
        OntSeqServiceClient.WorkspaceUri(
            18767, DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId).AbsoluteUri,
        "live workspace carries the exact Desktop-selected profile without a credential");
    AssertThrows<InvalidOperationException>(
        () => OntSeqServiceClient.WorkspaceUri(18767, "AML_LCWGS_UNKNOWN"),
        "workspace URI refuses an unknown profile handoff");
    AssertEqual(
        "http://127.0.0.1:18767/workspace?profile=AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25",
        OntSeqServiceClient.WorkspaceUri(
            18767,
            DesktopProfiles.Grch37UcscHg19Canonical25AdaptiveSamplingProfileId).AbsoluteUri,
        "live workspace carries the exact GRCh37 hg19 Adaptive Sampling profile");

    var canonicalOnlyConfig = new ServiceConfigResponse(
        "0.7.1", "/tmp/results", false, [],
        [DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId]);
    AssertEqual("True", canonicalOnlyConfig.SupportsProfile(
        DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId).ToString(),
        "service config resolves the exact requested Canonical-25 profile");
    AssertEqual("False", canonicalOnlyConfig.SupportsProfile(
        DesktopProfiles.DefaultProfileId).ToString(),
        "service config never substitutes its first or default profile");
    AssertEqual(
        DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId,
        canonicalOnlyConfig.RequireProfile(
            DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId).Profiles?.Single(),
        "service reuse accepts the selected profile when explicitly advertised");
    try
    {
        canonicalOnlyConfig.RequireProfile(DesktopProfiles.DefaultProfileId);
        throw new InvalidOperationException("service profile mismatch: expected refusal");
    }
    catch (InvalidOperationException error)
    {
        AssertEqual("True", error.Message.Contains(
            DesktopProfiles.DefaultProfileId, StringComparison.Ordinal).ToString(),
            "service profile mismatch names the exact requested profile");
    }

    var lcwgsOnly37Config = new ServiceConfigResponse(
        "0.7.1", "/tmp/results", false, [], [DesktopProfiles.Grch37LcwgsProfileId]);
    AssertEqual(
        DesktopProfiles.Grch37LcwgsProfileId,
        lcwgsOnly37Config.RequireProfile(DesktopProfiles.Grch37LcwgsProfileId).Profiles?.Single(),
        "GRCh37 lcWGS remains independently resolvable without an advertised panel profile");
    AssertThrows<InvalidOperationException>(
        () => lcwgsOnly37Config.RequireProfile(DesktopProfiles.Grch37AdaptiveSamplingProfileId),
        "GRCh37 Adaptive Sampling is not substituted when only lcWGS resolves");
    var adaptiveOnly37Config = new ServiceConfigResponse(
        "0.7.1", "/tmp/results", false, [], [DesktopProfiles.Grch37AdaptiveSamplingProfileId]);
    AssertEqual(
        DesktopProfiles.Grch37AdaptiveSamplingProfileId,
        adaptiveOnly37Config.RequireProfile(
            DesktopProfiles.Grch37AdaptiveSamplingProfileId).Profiles?.Single(),
        "GRCh37 Adaptive Sampling resolves only by its exact advertised identifier");

    const string expectedInstanceId = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    var launchExpectation = new ServiceLaunchExpectation(
        expectedInstanceId,
        "/opt/ontseq/resources",
        "/mnt/c/results",
        "/mnt/c/inputs");
    var matchingLaunchConfig = new ServiceConfigResponse(
        "0.7.1", "/mnt/c/results", false, [],
        [DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId],
        expectedInstanceId,
        "/opt/ontseq/resources",
        [new ServiceRootResponse("/mnt/c/inputs", @"C:\inputs")]);
    AssertEqual(
        expectedInstanceId,
        matchingLaunchConfig.RequireLaunch(launchExpectation).InstanceId,
        "service bootstrap accepts only the exact per-launch identity and filesystem scope");
    AssertThrows<ServiceLaunchScopeMismatchException>(
        () => matchingLaunchConfig.RequireLaunch(
            launchExpectation with { OutputDir = "/mnt/c/other-results" }),
        "service bootstrap rejects a mismatched output scope");
    AssertThrows<ServiceLaunchScopeMismatchException>(
        () => matchingLaunchConfig.RequireLaunch(
            launchExpectation with { AllowedRoot = "/mnt/c/other-inputs" }),
        "service bootstrap rejects a mismatched input scope");
    AssertThrows<ServiceInstanceMismatchException>(
        () => matchingLaunchConfig.RequireLaunch(
            launchExpectation with { InstanceId = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb" }),
        "service bootstrap rejects a foreign instance nonce");

    await using (var foreignService = new FakeOntSeqService(
        "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        launchExpectation.ResourceRoot,
        launchExpectation.OutputDir,
        launchExpectation.AllowedRoot))
    using (var foreignClient = new OntSeqServiceClient(foreignService.Port))
    {
        await AssertThrowsAsync<ServiceInstanceMismatchException>(
            () => foreignClient.BootstrapAsync(
                TimeSpan.FromSeconds(3), () => false, () => "", launchExpectation,
                CancellationToken.None),
            "Desktop never bootstraps an already-listening foreign service");
        AssertThrows<InvalidOperationException>(
            () => foreignClient.GetConfigAsync(CancellationToken.None).GetAwaiter().GetResult(),
            "rejected foreign client remains unbootstrapped and cannot be reused");
    }

    await using (var wrongScopeService = new FakeOntSeqService(
        expectedInstanceId,
        launchExpectation.ResourceRoot,
        "/mnt/c/foreign-results",
        launchExpectation.AllowedRoot))
    using (var wrongScopeClient = new OntSeqServiceClient(wrongScopeService.Port))
    {
        await AssertThrowsAsync<ServiceLaunchScopeMismatchException>(
            () => wrongScopeClient.BootstrapAsync(
                TimeSpan.FromSeconds(3), () => false, () => "", launchExpectation,
                CancellationToken.None),
            "matching nonce cannot hide a mismatched launch scope");
        AssertThrows<InvalidOperationException>(
            () => wrongScopeClient.GetConfigAsync(CancellationToken.None).GetAwaiter().GetResult(),
            "scope-rejected client remains unbootstrapped and cannot be reused");
    }
    AssertEqual("True", WslServiceLauncher.IsServiceLaunchCollision(
            new ServiceInstanceMismatchException(
                expectedInstanceId, "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb")).ToString(),
        "foreign listener triggers isolated-port retry");
    AssertEqual("True", WslServiceLauncher.IsServiceLaunchCollision(
            new InvalidOperationException("Backend exited"), "OSError: [Errno 98] Address already in use").ToString(),
        "backend bind collision triggers isolated-port retry");
    AssertEqual("False", WslServiceLauncher.IsServiceLaunchCollision(
            new InvalidOperationException("Reference bundle missing"), "profile incomplete").ToString(),
        "resource failures do not masquerade as port collisions");

    using (var occupiedListener = new System.Net.Sockets.TcpListener(
               System.Net.IPAddress.Loopback, 0))
    {
        occupiedListener.Server.ExclusiveAddressUse = true;
        occupiedListener.Start();
        var occupiedPort = ((System.Net.IPEndPoint)occupiedListener.LocalEndpoint).Port;
        var isolatedPort = WslServiceLauncher.SelectLoopbackPort(occupiedPort);
        AssertEqual("True", isolatedPort.UsedFallback.ToString(),
            "occupied configured port selects an isolated fallback");
        AssertEqual("False", (isolatedPort.Port == occupiedPort).ToString(),
            "new Desktop never attaches to an existing listener");
        AssertEqual("True", (isolatedPort.Port is >= 1 and <= 65535).ToString(),
            "isolated fallback is a valid TCP port");
        var launchSettings = new DesktopSettings
        {
            ResourceRootWsl = "/opt/ontseq/resources",
            OutputDirectoryWindows = Path.Combine(root, "isolated-port-results"),
            Port = occupiedPort
        };
        var serviceInstanceId = WslServiceLauncher.CreateServiceInstanceId();
        var launchArguments = WslServiceLauncher.ProfileServiceArguments(
            launchSettings, root, isolatedPort.Port, serviceInstanceId).ToArray();
        var portArgument = Array.IndexOf(launchArguments, "--port");
        AssertEqual(isolatedPort.Port.ToString(), launchArguments[portArgument + 1],
            "WSL service receives the isolated port rather than the occupied setting");
        var instanceArgument = Array.IndexOf(launchArguments, "--instance-id");
        AssertEqual(serviceInstanceId, launchArguments[instanceArgument + 1],
            "WSL service receives the Desktop's unique per-launch identity");
    }
    AssertThrows<ArgumentOutOfRangeException>(
        () => WslServiceLauncher.SelectLoopbackPort(0),
        "invalid configured service port refusal");
    AssertEqual(root, WslServiceLauncher.WorkspaceAllowedRootWindows(bam, Path.Combine(root, "results")),
        "workspace grants only selected BAM parent");
    AssertEqual(Path.Combine(root, "results"), WslServiceLauncher.WorkspaceAllowedRootWindows("", Path.Combine(root, "results")),
        "workspace without BAM grants only output directory");
    AssertEqual("True", WslServiceLauncher.IsWithinAllowedRoot(Path.Combine(root, "child"), root).ToString(),
        "existing service can safely reuse a nested input folder");
    AssertEqual("False", WslServiceLauncher.IsWithinAllowedRoot(root + "-other", root).ToString(),
        "existing service never expands scope via path-prefix collision");

    var runtimeFolder = Path.Combine(root, "runtime");
    Directory.CreateDirectory(runtimeFolder);
    var baseArchive = Path.Combine(runtimeFolder, "ontseq-linux-runtime.tar.gz");
    var coreWheel = Path.Combine(runtimeFolder, "ontseq_platform-0.7.1-py3-none-any.whl");
    var checksumsPath = Path.Combine(runtimeFolder, "SHA256SUMS");
    File.WriteAllText(baseArchive, "synthetic base archive: not executable");
    File.WriteAllText(coreWheel, "synthetic Core wheel: not executable");
    var archiveChecksum = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(baseArchive)));
    var wheelChecksum = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(coreWheel)));
    File.WriteAllLines(checksumsPath, [
        archiveChecksum + "  " + Path.GetFileName(baseArchive),
        wheelChecksum + "  " + Path.GetFileName(coreWheel)
    ]);
    var package = await RuntimePackage.VerifyAsync(baseArchive, "0.7.1", CancellationToken.None);
    AssertEqual(coreWheel, package.WheelPath, "runtime requires matching Core wheel sidecar");
    var installCommand = WslServiceLauncher.RuntimeInstallCoreCommand(
        package, "/home/synthetic/.local/share/ontseq/runtime-v0.7.1-synthetic");
    AssertEqual("False", installCommand.Contains("rm ", StringComparison.Ordinal).ToString(),
        "runtime installation never removes an existing prefix");
    AssertEqual("True", installCommand.StartsWith("test ! -e ", StringComparison.Ordinal).ToString(),
        "runtime installation requires a new prefix");
    AssertEqual("True", (installCommand.IndexOf("/conda-unpack", StringComparison.Ordinal) <
        installCommand.IndexOf("-m pip", StringComparison.Ordinal)).ToString(),
        "base relocation precedes wheel installation");
    AssertEqual("True", installCommand.Contains("--no-deps --no-index --force-reinstall", StringComparison.Ordinal).ToString(),
        "Core update is strictly offline without changing dependency packages");
    AssertEqual("True", installCommand.Contains("assert value ==", StringComparison.Ordinal).ToString(),
        "runtime Core version must match before settings can be saved");
    AssertEqual("False", installCommand.Contains("$(", StringComparison.Ordinal).ToString(),
        "version verification cannot be evaluated early by an outer WSL shell");
    AssertThrows<InvalidDataException>(() => WslServiceLauncher.RuntimeInstallCoreCommand(package, "/home/synthetic"),
        "runtime installation rejects broad or old target paths");
    AssertThrows<InvalidDataException>(() => WslServiceLauncher.RuntimeInstallCoreCommand(
        package, "/home/synthetic/.local/share/ontseq/runtime-v0.5.3-existing"),
        "runtime installation cannot overwrite old release prefix");
    File.AppendAllText(coreWheel, "tampered");
    AssertThrows<InvalidDataException>(() => RuntimePackage.VerifyAsync(baseArchive, "0.7.1", CancellationToken.None)
        .GetAwaiter().GetResult(), "tampered Core wheel blocks installation");
    File.WriteAllText(coreWheel, "synthetic Core wheel: not executable");
    File.AppendAllText(baseArchive, "tampered");
    AssertThrows<InvalidDataException>(() => RuntimePackage.VerifyAsync(baseArchive, "0.7.1", CancellationToken.None)
        .GetAwaiter().GetResult(), "tampered base archive blocks installation");
    File.WriteAllText(baseArchive, "synthetic base archive: not executable");
    File.WriteAllText(checksumsPath, archiveChecksum + "  " + Path.GetFileName(baseArchive));
    AssertThrows<InvalidDataException>(() => RuntimePackage.VerifyAsync(baseArchive, "0.7.1", CancellationToken.None)
        .GetAwaiter().GetResult(), "missing Core checksum blocks installation");

    var originalSettingsOverride = Environment.GetEnvironmentVariable("ONTSEQ_DESKTOP_SETTINGS");
    try
    {
        var isolatedSettingsPath = Path.Combine(root, "isolated", "desktop.settings.json");
        Environment.SetEnvironmentVariable("ONTSEQ_DESKTOP_SETTINGS", isolatedSettingsPath);
        AssertEqual(isolatedSettingsPath, DesktopSettings.UserSettingsPath,
            "isolated Desktop shows and writes only its explicit settings path");
        AssertEqual(DesktopProfiles.DefaultProfileId, DesktopSettings.Load().DefaultProfile,
            "missing isolated settings use defaults without falling back to existing installation");
        new DesktopSettings {
            DefaultProfile = DesktopProfiles.Grch37LcwgsProfileId,
            ModkitExecutableWsl = separateModkitPath
        }.SaveUserSettings();
        AssertEqual("True", File.Exists(isolatedSettingsPath).ToString(), "isolated settings save destination");
        AssertEqual(DesktopProfiles.Grch37LcwgsProfileId, DesktopSettings.Load().DefaultProfile,
            "isolated Desktop loads its own native37 profile");
        AssertEqual(separateModkitPath, DesktopSettings.Load().ModkitExecutableWsl,
            "optional modkit executable survives isolated settings serialization and load");
        foreach (var invalid in new[] { "relative.json", @"\\server\share\settings.json", Path.Combine(root, "settings.cmd") })
        {
            Environment.SetEnvironmentVariable("ONTSEQ_DESKTOP_SETTINGS", invalid);
            AssertThrows<InvalidDataException>(() => DesktopSettings.Load(),
                "invalid settings override fails closed instead of using global settings");
        }
    }
    finally
    {
        Environment.SetEnvironmentVariable("ONTSEQ_DESKTOP_SETTINGS", originalSettingsOverride);
    }

    const string legacySettingsJson = """
        {
          "wslDistribution": "Ubuntu",
          "referenceLocksWsl": {"GRCh37": "/legacy/grch37.lock.json"},
          "adaptiveTargetBedWsl": "/legacy/roi.bed",
          "adaptiveTargetBedVersion": "legacy-v1"
        }
        """;
    var legacySettings = JsonSerializer.Deserialize<DesktopSettings>(
        legacySettingsJson, JsonDefaults.Options)
        ?? throw new InvalidOperationException("legacy settings deserialization returned null");
    legacySettings.ApplyProfileDefaults();
    AssertEqual(
        "/legacy/grch37.lock.json",
        legacySettings.ReferenceLocksWsl["GRCh37"],
        "legacy explicit reference remains readable");
    AssertEqual("/legacy/roi.bed", legacySettings.AdaptiveTargetBedWsl, "legacy BED remains readable");
    AssertEqual(
        profileDefaults.ResourceRootWsl,
        legacySettings.ResourceRootWsl,
        "legacy settings gain user-writable resource root");
    AssertEqual(
        "AML_LCWGS_GRCh38",
        legacySettings.DefaultProfile,
        "legacy settings gain GRCh38 default profile");

    var requestJson = JsonSerializer.Serialize(
        new RunStartRequest(
            @"C:\Lab\sample.bam",
            "SAMPLE_001",
            null,
            "AML_AS_111_GRCh38",
            "GRCh38",
            "adaptive_sampling"),
        JsonDefaults.Options);
    using (var requestDocument = JsonDocument.Parse(requestJson))
    {
        AssertEqual(
            "AML_AS_111_GRCh38",
            requestDocument.RootElement.GetProperty("profile").GetString(),
            "profile API field");
        AssertEqual(
            "GRCh38",
            requestDocument.RootElement.GetProperty("genome_build").GetString(),
            "GRCh38 compatibility API field");
        AssertEqual(
            "False",
            requestDocument.RootElement.TryGetProperty("run_id", out _).ToString(),
            "run ID is omitted so Core derives sample plus UTC timestamp");
    }

    var grch37AdaptiveRequestJson = JsonSerializer.Serialize(
        new RunStartRequest(
            @"C:\Lab\sample.bam",
            "SAMPLE_037_AS",
            null,
            DesktopProfiles.Grch37UcscHg19Canonical25AdaptiveSamplingProfileId,
            "GRCh37",
            "adaptive_sampling"),
        JsonDefaults.Options);
    using (var requestDocument = JsonDocument.Parse(grch37AdaptiveRequestJson))
    {
        AssertEqual(
            DesktopProfiles.Grch37UcscHg19Canonical25AdaptiveSamplingProfileId,
            requestDocument.RootElement.GetProperty("profile").GetString(),
            "GRCh37 Adaptive Sampling profile API field");
        AssertEqual("GRCh37",
            requestDocument.RootElement.GetProperty("genome_build").GetString(),
            "GRCh37 Adaptive Sampling build API field");
        AssertEqual("adaptive_sampling",
            requestDocument.RootElement.GetProperty("assay").GetString(),
            "GRCh37 Adaptive Sampling assay API field");
    }

    const string upperHash = "E518E7131D51ABED37A7AED5DB6A031B753ADF424E867B52352B44CA4A6E7B4B";
    AssertEqual(
        "GRCh38_LOCAL_e518e7131d51abed",
        WslServiceLauncher.ReferenceIdFor("GRCh38", upperHash),
        "content-addressed reference ID");
    AssertEqual(
        "GRCh38.e518e7131d51abed.reference-lock.json",
        WslServiceLauncher.ReferenceLockFileNameFor("GRCh38", upperHash),
        "content-addressed reference filename");
    AssertThrows<ArgumentException>(
        () => WslServiceLauncher.ReferenceIdFor("GRCh38", "not-a-sha256"),
        "invalid reference fingerprint refusal");

    var existingLock = Path.Combine(root, "existing.reference-lock.json");
    var rejectedLock = Path.Combine(root, "rejected.reference-lock.tmp.json");
    File.WriteAllText(existingLock, "existing-lock");
    File.WriteAllText(rejectedLock, "{\"source_fai_sha256\":\"" + new string('b', 64) + "\"}");
    AssertThrows<InvalidDataException>(
        () => WslServiceLauncher.PublishReferenceLockFile(
            rejectedLock,
            existingLock,
            new string('a', 64)),
        "changed FAI refusal");
    AssertEqual("existing-lock", File.ReadAllText(existingLock), "active lock preservation");

    var adaptiveSettings = new DesktopSettings
    {
        AdaptiveTargetBedWsl = "/mnt/c/ONTSeq/resources/adaptive_sampling/analysis_roi.test.bed",
        AdaptiveTargetBedVersion = "panel.bed@sha256:" + new string('a', 64)
    };
    AssertEqual("True", adaptiveSettings.HasAdaptiveTargetBedConfiguration.ToString(), "adaptive BED configured state");
    adaptiveSettings.ClearAdaptiveTargetBed();
    AssertEqual(null, adaptiveSettings.AdaptiveTargetBedWsl, "adaptive BED path cleared");
    AssertEqual(null, adaptiveSettings.AdaptiveTargetBedVersion, "adaptive BED version cleared");
    AssertEqual("False", adaptiveSettings.HasAdaptiveTargetBedConfiguration.ToString(), "adaptive BED cleared state");

    var partialAdaptiveSettings = new DesktopSettings
    {
        AdaptiveTargetBedVersion = "incomplete-bed-version"
    };
    AssertEqual("True", partialAdaptiveSettings.HasAdaptiveTargetBedConfiguration.ToString(), "partial adaptive BED remains removable");
    partialAdaptiveSettings.ClearAdaptiveTargetBed();
    AssertEqual("False", partialAdaptiveSettings.HasAdaptiveTargetBedConfiguration.ToString(), "partial adaptive BED clear state");

    Console.WriteLine("Desktop profile, resource bridge, compatibility, path and BAM index tests passed.");
}
finally
{
    Directory.Delete(root, recursive: true);
}

static void AssertEqual(string? expected, string? actual, string scenario)
{
    if (!string.Equals(expected, actual, StringComparison.Ordinal))
    {
        throw new InvalidOperationException(
            $"{scenario}: expected '{expected ?? "<null>"}', got '{actual ?? "<null>"}'.");
    }
}

static void AssertThrows<TException>(Action action, string scenario) where TException : Exception
{
    try
    {
        action();
    }
    catch (TException)
    {
        return;
    }
    throw new InvalidOperationException($"{scenario}: expected {typeof(TException).Name}.");
}

static async Task AssertThrowsAsync<TException>(Func<Task> action, string scenario)
    where TException : Exception
{
    try
    {
        await action();
    }
    catch (TException)
    {
        return;
    }
    throw new InvalidOperationException($"{scenario}: expected {typeof(TException).Name}.");
}

static void AssertSequenceEqual(
    IReadOnlyList<string> expected,
    IReadOnlyList<string> actual,
    string scenario)
{
    if (!expected.SequenceEqual(actual, StringComparer.Ordinal))
    {
        throw new InvalidOperationException(
            $"{scenario}: expected '{string.Join(" ", expected)}', got '{string.Join(" ", actual)}'.");
    }
}

static async Task VerifyMethylationHttpContractAsync()
{
    const string bamPath = "/synthetic/BAM mit Umlaut ö.bam";
    const string scanId = "0123456789abcdef0123456789abcdef";
    var probe = new MethylationProbeResponse(bamPath, "detected", "Synthetic 5mC witness", 42000, false,
        true, ReasonCode: "detected_5mc", ElapsedSeconds: 4.2, ScanMode: "thorough", Reader: "pysam");
    var requests = new System.Collections.Concurrent.ConcurrentQueue<FakeServiceRequest>();
    var scanStarts = 0;
    var progressReads = 0;
    await using var server = new FakeOntSeqService(new string('b', 32), "/synthetic/resources", "/synthetic/output", "/synthetic",
        request =>
        {
            requests.Enqueue(request);
            if (request.Path == "/api/methylation/probe") return new(200, JsonSerializer.Serialize(probe, JsonDefaults.Options));
            if (request.Path == "/api/methylation/scans")
            {
                if (++scanStarts > 1) return new(409, "{\"error\":\"Synthetic busy\",\"reason_code\":\"busy\"}");
                return new(202, JsonSerializer.Serialize(new MethylationScanSnapshot(scanId, bamPath, "running", 0, 0, null), JsonDefaults.Options));
            }
            if (request.Path == $"/api/methylation/scans/{scanId}/cancel")
                return new(200, JsonSerializer.Serialize(new MethylationScanSnapshot(scanId, bamPath, "cancelled", 5000, 1.2,
                    probe with { Status = "unknown", ReasonCode = "cancelled", CheckedReads = 5000, ElapsedSeconds = 1.2 }), JsonDefaults.Options));
            if (request.Path == $"/api/methylation/scans/{scanId}")
                return new(200, JsonSerializer.Serialize(new MethylationScanSnapshot(scanId,
                    ++progressReads == 1 ? bamPath : "/synthetic/wrong.bam", "running", 5000, 1.2, null), JsonDefaults.Options));
            return null;
        });
    using var client = new OntSeqServiceClient(server.Port);
    await client.BootstrapAsync(TimeSpan.FromSeconds(5), null, null,
        new ServiceLaunchExpectation(new string('b', 32), "/synthetic/resources", "/synthetic/output", "/synthetic"), CancellationToken.None);
    var fresh = await client.ProbeMethylationAsync(bamPath, CancellationToken.None, forceRefresh: true);
    AssertEqual("detected", fresh.Status, "fresh prestart probe preserves a server-confirmed thorough witness");
    var started = await client.StartMethylationScanAsync(bamPath, CancellationToken.None);
    AssertEqual("running", started.State, "202 scan creation returns immediately with typed state");
    var progress = await client.GetMethylationScanAsync(started.ScanId, bamPath, CancellationToken.None);
    AssertEqual("5000", progress.CheckedReads.ToString(), "asynchronous progress comes from the server");
    var cancelled = await client.CancelMethylationScanAsync(started.ScanId, bamPath, CancellationToken.None);
    AssertEqual("cancelled", cancelled.State, "cancel uses the concrete scan identity");
    AssertEqual("unknown", cancelled.Result?.Status, "cancel never implies absent methylation");
    try
    {
        await client.StartMethylationScanAsync(bamPath, CancellationToken.None);
        throw new InvalidOperationException("Expected synthetic busy response");
    }
    catch (MethylationServiceException error)
    {
        AssertEqual("busy", error.ReasonCode, "HTTP errors retain the structured reason for German presentation");
    }
    await AssertThrowsAsync<InvalidDataException>(
        () => client.GetMethylationScanAsync(scanId, bamPath, CancellationToken.None),
        "client refuses progress from a different selected BAM");
    var calls = requests.Where(request => request.Path.StartsWith("/api/methylation/", StringComparison.Ordinal)).ToArray();
    AssertEqual("POST", calls[0].Method, "prestart probe uses POST");
    using (var json = JsonDocument.Parse(calls[0].Body))
    {
        AssertEqual("True", json.RootElement.GetProperty("force_refresh").GetBoolean().ToString(),
            "actual prestart HTTP request bypasses stale cache");
        AssertEqual(bamPath, json.RootElement.GetProperty("bam_path").GetString(), "UTF-8 BAM path is transported opaquely");
    }
    AssertEqual("True", calls.All(request => request.Headers.Contains("X-ONTSeq-Token: synthetic-service-token", StringComparison.OrdinalIgnoreCase)).ToString(),
        "scan endpoints preserve authenticated local-service access");
    AssertEqual("{}", calls.Single(request => request.Path.EndsWith("/cancel", StringComparison.Ordinal)).Body,
        "cancel sends an explicit empty JSON object");
    AssertEqual("False", calls.Any(request => request.Path == "/api/runs").ToString(),
        "probing and cancellation do not start an analysis");
}

static void VerifyMethylationLayout()
{
    Exception? failure = null;
    var thread = new Thread(() =>
    {
        try
        {
            // Layout only: never Show(), raise Loaded, load settings, or inspect input files.
            var window = new MainWindow();
            var content = (System.Windows.FrameworkElement)window.Content;
            var explanation = (System.Windows.Controls.TextBlock)window.FindName("MethylationExplanationText");
            explanation.Text = string.Join(" ", Enumerable.Repeat(
                "Die BAM wurde noch nicht vollständig geprüft; es ist keine sichere Aussage über vorhandene 5mC-Tags möglich.", 8));
            var details = (System.Windows.Controls.TextBlock)window.FindName("MethylationDetailsText");
            details.Text = string.Join("\n", Enumerable.Repeat("Synthetic read_error: zusätzliche verständliche Prüfdetails.", 12));
            ((System.Windows.Controls.Expander)details.Parent).IsExpanded = true;
            content.Measure(new System.Windows.Size(944, 720));
            content.Arrange(new System.Windows.Rect(0, 0, 944, 720));
            content.UpdateLayout();
            var scroll = (System.Windows.Controls.ScrollViewer)window.FindName("InputsScroll");
            AssertEqual("True", (scroll.ScrollableHeight > 0).ToString(), "long methylation details scroll at minimum window size");
            var analysis = (System.Windows.FrameworkElement)window.FindName("AnalysisProgressPanel");
            AssertEqual("True", (analysis.ActualHeight >= 100).ToString(), "analysis progress retains usable space beside long details");
            scroll.ScrollToBottom();
            content.UpdateLayout();
            var start = (System.Windows.FrameworkElement)window.FindName("StartButton");
            var point = start.TransformToAncestor(scroll).Transform(new System.Windows.Point(0, 0));
            AssertEqual("True", (point.Y >= -1 && point.Y + start.ActualHeight <= scroll.ActualHeight + 1).ToString(),
                "analysis start remains reachable through form scrolling");
        }
        catch (Exception error) { failure = error; }
    });
    thread.SetApartmentState(ApartmentState.STA);
    thread.Start();
    thread.Join();
    if (failure is not null) throw new InvalidOperationException("WPF methylation layout check failed", failure);
}

sealed record FakeServiceRequest(string Method, string Path, string Body, string Headers);
sealed record FakeServiceReply(int StatusCode, string Body);

sealed class FakeOntSeqService : IAsyncDisposable
{
    private readonly TcpListener _listener = new(IPAddress.Loopback, 0);
    private readonly CancellationTokenSource _cancellation = new();
    private readonly Task _serverTask;
    private readonly string _configJson;
    private readonly Func<FakeServiceRequest, FakeServiceReply?>? _handler;

    public FakeOntSeqService(
        string instanceId,
        string resourceRoot,
        string outputDir,
        string allowedRoot,
        Func<FakeServiceRequest, FakeServiceReply?>? handler = null)
    {
        _handler = handler;
        _configJson = JsonSerializer.Serialize(new
        {
            version = "0.7.1",
            output_dir = outputDir,
            busy = false,
            not_wired = Array.Empty<string>(),
            profiles = new[] { DesktopProfiles.Grch37UcscHg19Canonical25LcwgsProfileId },
            instance_id = instanceId,
            resource_root = resourceRoot,
            roots = new[] { new { posix = allowedRoot, display = allowedRoot } }
        });
        _listener.Server.ExclusiveAddressUse = true;
        _listener.Start();
        Port = ((IPEndPoint)_listener.LocalEndpoint).Port;
        _serverTask = ServeAsync(_cancellation.Token);
    }

    public int Port { get; }

    private async Task ServeAsync(CancellationToken cancellationToken)
    {
        try
        {
            while (!cancellationToken.IsCancellationRequested)
            {
                using var client = await _listener.AcceptTcpClientAsync(cancellationToken);
                await RespondAsync(client, cancellationToken);
            }
        }
        catch (OperationCanceledException)
        {
        }
        catch (SocketException) when (cancellationToken.IsCancellationRequested)
        {
        }
    }

    private async Task RespondAsync(TcpClient client, CancellationToken cancellationToken)
    {
        using var stream = client.GetStream();
        using var requestBytes = new MemoryStream();
        var buffer = new byte[2048];
        var headerEnd = -1;
        while (headerEnd < 0)
        {
            var count = await stream.ReadAsync(buffer, cancellationToken);
            if (count == 0) break;
            requestBytes.Write(buffer, 0, count);
            headerEnd = Encoding.ASCII.GetString(requestBytes.ToArray()).IndexOf("\r\n\r\n", StringComparison.Ordinal);
            if (requestBytes.Length > 16384) throw new InvalidDataException("Synthetic HTTP request too large.");
        }
        if (headerEnd < 0) return;
        var headers = Encoding.ASCII.GetString(requestBytes.ToArray(), 0, headerEnd);
        var lengthHeader = headers.Split("\r\n").FirstOrDefault(line => line.StartsWith("Content-Length:", StringComparison.OrdinalIgnoreCase));
        var contentLength = lengthHeader is null ? 0 : int.Parse(lengthHeader.Split(':')[1].Trim());
        while (requestBytes.Length < headerEnd + 4 + contentLength)
        {
            var count = await stream.ReadAsync(buffer, cancellationToken);
            if (count == 0) break;
            requestBytes.Write(buffer, 0, count);
        }
        var firstLine = headers.Split("\r\n")[0].Split(' ');
        var target = firstLine[1];
        string requestBody;
        if (headers.Contains("Transfer-Encoding: chunked", StringComparison.OrdinalIgnoreCase))
        {
            using var decoded = new MemoryStream();
            var offset = headerEnd + 4;
            while (true)
            {
                var lineEnd = Encoding.ASCII.GetString(requestBytes.ToArray()).IndexOf("\r\n", offset, StringComparison.Ordinal);
                while (lineEnd < 0)
                {
                    await ReadMoreAsync();
                    lineEnd = Encoding.ASCII.GetString(requestBytes.ToArray()).IndexOf("\r\n", offset, StringComparison.Ordinal);
                }
                var length = int.Parse(Encoding.ASCII.GetString(requestBytes.ToArray(), offset, lineEnd - offset),
                    System.Globalization.NumberStyles.HexNumber);
                offset = lineEnd + 2;
                if (length == 0) break;
                while (requestBytes.Length < offset + length + 2) await ReadMoreAsync();
                decoded.Write(requestBytes.ToArray(), offset, length);
                offset += length + 2;
            }
            requestBody = Encoding.UTF8.GetString(decoded.ToArray());
        }
        else requestBody = Encoding.UTF8.GetString(requestBytes.ToArray(), headerEnd + 4, contentLength);
        async Task ReadMoreAsync()
        {
            var count = await stream.ReadAsync(buffer, cancellationToken);
            if (count == 0) throw new EndOfStreamException("Incomplete synthetic request body");
            requestBytes.Write(buffer, 0, count);
            if (requestBytes.Length > 16384) throw new InvalidDataException("Synthetic HTTP request too large.");
        }
        var custom = _handler?.Invoke(new FakeServiceRequest(firstLine[0], target, requestBody, headers));
        var body = custom?.Body ?? (target == "/api/config"
            ? _configJson
            : "<!doctype html><script>const TOKEN = \"synthetic-service-token\";</script>");
        var contentType = target.StartsWith("/api/", StringComparison.Ordinal) ? "application/json" : "text/html";
        var payload = Encoding.UTF8.GetBytes(body);
        var header = Encoding.ASCII.GetBytes(
            $"HTTP/1.1 {custom?.StatusCode ?? 200} Synthetic\r\nContent-Type: {contentType}\r\nContent-Length: {payload.Length}\r\n" +
            "Connection: close\r\nCache-Control: no-store\r\n\r\n");
        await stream.WriteAsync(header, cancellationToken);
        await stream.WriteAsync(payload, cancellationToken);
    }

    public async ValueTask DisposeAsync()
    {
        _cancellation.Cancel();
        _listener.Stop();
        await _serverTask;
        _cancellation.Dispose();
    }
}
