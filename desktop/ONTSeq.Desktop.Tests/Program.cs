using ONTSeq.Desktop;

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

    Console.WriteLine("Desktop path and BAM index tests passed.");
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
