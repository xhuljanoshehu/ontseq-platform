namespace ONTSeq.Desktop;

public static class BamIndexLocator
{
    public static string? Find(string bamPath)
    {
        if (string.IsNullOrWhiteSpace(bamPath)) return null;

        var preferred = bamPath + ".bai";
        if (File.Exists(preferred)) return preferred;

        var alternative = Path.ChangeExtension(bamPath, ".bai");
        return File.Exists(alternative) ? alternative : null;
    }
}
