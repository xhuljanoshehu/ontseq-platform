namespace ONTSeq.Desktop;

/// <summary>A preparation failure attributed to one dependency, before any analysis starts.</summary>
public sealed class WslPrerequisiteException : InvalidOperationException
{
    public string ReasonCode { get; }
    public string? WslPath { get; }
    public string? WindowsPath { get; }

    public WslPrerequisiteException(
        string reasonCode,
        string message,
        string? wslPath = null,
        string? windowsPath = null,
        Exception? innerException = null) : base(message, innerException)
    {
        ReasonCode = reasonCode;
        WslPath = wslPath;
        WindowsPath = windowsPath;
    }
}

internal sealed record WslPrerequisiteCheck(
    string ReasonCode,
    string WslPath,
    string? WindowsPath,
    string FailureMessage,
    string TestCommand)
{
    public WslPrerequisiteException Failure() =>
        new(ReasonCode, FailureMessage, WslPath, WindowsPath);
}
