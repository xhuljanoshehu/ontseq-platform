using ONTSeq.Desktop;
using System.Text.Json;

var root = Path.Combine(Path.GetTempPath(), "ONTSeq.Desktop.Tests", Guid.NewGuid().ToString("N"));
Directory.CreateDirectory(root);

try
{
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
        "~/.local/share/ontseq/resources",
        profileDefaults.ResourceRootWsl,
        "user-writable default WSL resource root");
    AssertEqual("AML_LCWGS_GRCh38", profileDefaults.DefaultProfile, "default analysis profile");
    AssertEqual("4", DesktopProfiles.Supported.Count.ToString(), "exact supported profile count");
    AssertEqual(
        "True",
        DesktopProfiles.Supported.All(profile => profile.GenomeBuild == "GRCh38").ToString(),
        "all published desktop profiles are GRCh38-only");
    AssertEqual(
        "AML_AS_111_GRCh38",
        DesktopProfiles.Require("AML_AS_111_GRCh38").ProfileId,
        "adaptive sampling profile identity");
    AssertEqual(
        "Canonical-25 (chr1–22, chrX, chrY, chrM)",
        DesktopProfiles.Require("AML_LCWGS_GRCh38_CANONICAL25").DictionaryLabel,
        "Canonical-25 profile exposes its exact BAM dictionary contract");
    AssertThrows<InvalidOperationException>(
        () => DesktopProfiles.Require("AML_AS_111_GRCh37"),
        "GRCh37 profile refusal");
    AssertSequenceEqual(
        [
            "GRCh38_GENCODE50_MANE1.5_v1",
            "HEMATOLOGY_v1",
            "AML_AS_111_GRCh38_v1"
        ],
        WslServiceLauncher.ManagedGrch38ResourceBundleIds,
        "Desktop repair owns the complete GRCh38 profile resource family");

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

    AssertSequenceEqual(
        ["references", "status", "--resource-root", "/opt/ontseq"],
        WslServiceLauncher.ResourceManagementArguments("status", "/opt/ontseq/"),
        "bundle status command bridge");
    AssertSequenceEqual(
        [
            "references", "status", "--resource-root",
            "~/.local/share/ontseq/resources"
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
        "~/.local/share/ontseq/resources",
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
