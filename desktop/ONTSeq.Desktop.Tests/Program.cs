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

    Console.WriteLine("BAM index detection tests passed.");
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
