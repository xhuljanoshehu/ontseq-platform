using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using Microsoft.Win32;

namespace ONTSeq.Desktop;

public partial class MainWindow : Window
{
    private static readonly Regex SafeId = new("^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$", RegexOptions.Compiled);

    private DesktopSettings _settings = new();
    private WslServiceLauncher? _launcher;
    private OntSeqServiceClient? _client;
    private CancellationTokenSource? _analysisCts;
    private string? _currentRunId;
    private string? _currentSampleId;

    public ObservableCollection<StageDisplay> StageItems { get; } = [];

    public MainWindow()
    {
        InitializeComponent();
        DataContext = this;
        SetPlaceholders();
    }

    private void Window_Loaded(object sender, RoutedEventArgs e)
    {
        ReloadSettingsState();
    }

    private void ReloadSettingsState()
    {
        try
        {
            _settings = DesktopSettings.Load();
            var configured = new List<string>();
            if (_settings.TryReferenceLockFor("GRCh38", out _)) configured.Add("GRCh38");
            if (_settings.TryReferenceLockFor("GRCh37", out _)) configured.Add("GRCh37");

            BackendStateText.Text = configured.Count == 0
                ? "Einrichtung erforderlich"
                : $"WSL: {_settings.WslDistribution} · Ref: {string.Join(", ", configured)}";

            if (configured.Count == 0)
                DetailText.Text = "Vor der ersten Analyse bitte 'System einrichten' öffnen: Linux-Runtime installieren, passende Referenz konfigurieren und Selbsttest ausführen.";
            StartButton.IsEnabled = true;
        }
        catch (Exception error)
        {
            BackendStateText.Text = "Konfigurationsfehler";
            DetailText.Text = error.Message;
            StartButton.IsEnabled = false;
        }
    }

    private void Setup_Click(object sender, RoutedEventArgs e)
    {
        var setup = new SetupWindow(_settings) { Owner = this };
        setup.ShowDialog();
        ReloadSettingsState();
    }

    private void Browse_Click(object sender, RoutedEventArgs e)
    {
        var dialog = new OpenFileDialog
        {
            Title = "Oxford-Nanopore BAM auswählen",
            Filter = "BAM-Dateien (*.bam)|*.bam",
            CheckFileExists = true,
            Multiselect = false
        };
        if (dialog.ShowDialog(this) != true) return;
        BamPathTextBox.Text = dialog.FileName;
        var sample = Path.GetFileNameWithoutExtension(dialog.FileName);
        SampleIdTextBox.Text = SanitizeSuggestedId(sample);
        var index = BamIndexLocator.Find(dialog.FileName);
        DetailText.Text = index is not null
            ? $"BAM und BAM-Index gefunden ({Path.GetFileName(index)})."
            : "Hinweis: Erwartet wird <sample>.bam.bai oder <sample>.bai neben der BAM-Datei.";
    }

    private async void Start_Click(object sender, RoutedEventArgs e)
    {
        if (!TryValidateForm(out var bam, out var sampleId, out var genomeBuild, out var assay)) return;

        if (!_settings.TryReferenceLockFor(genomeBuild, out var referenceLock))
        {
            RunStateText.Text = "EINRICHTUNG";
            DetailText.Text = $"Für {genomeBuild} fehlt der Reference-Lock. Die Analyse wurde noch nicht gestartet.";
            var openSetup = MessageBox.Show(
                this,
                $"Für {genomeBuild} ist noch keine passende Referenz eingerichtet. Jetzt 'System einrichten' öffnen?",
                "ONTSeq Einrichtung erforderlich",
                MessageBoxButton.YesNo,
                MessageBoxImage.Information);
            if (openSetup != MessageBoxResult.Yes) return;

            var setup = new SetupWindow(_settings) { Owner = this };
            setup.ShowDialog();
            ReloadSettingsState();
            if (!_settings.TryReferenceLockFor(genomeBuild, out referenceLock))
            {
                RunStateText.Text = "EINRICHTUNG";
                DetailText.Text = $"{genomeBuild}-Referenz ist weiterhin nicht konfiguriert; es wurde keine Analyse gestartet.";
                return;
            }
        }

        StartButton.IsEnabled = false;
        BrowseButton.IsEnabled = false;
        SetupButton.IsEnabled = false;
        OpenReportButton.IsEnabled = false;
        OpenExcelButton.IsEnabled = false;
        OpenFolderButton.IsEnabled = false;
        RunProgress.Value = 0;
        RunProgress.IsIndeterminate = true;
        SetPlaceholders();

        _analysisCts?.Cancel();
        _analysisCts?.Dispose();
        _analysisCts = new CancellationTokenSource();
        var cancellationToken = _analysisCts.Token;

        try
        {
            if (_client is not null) { _client.Dispose(); _client = null; }
            if (_launcher is not null) await _launcher.DisposeAsync();
            _launcher = new WslServiceLauncher();

            var allowedRoot = Path.GetDirectoryName(bam)
                              ?? throw new InvalidOperationException("BAM-Verzeichnis konnte nicht bestimmt werden.");

            BackendStateText.Text = "WSL, Backend und Referenz werden geprüft…";
            DetailText.Text = "Vorprüfung: WSL2, ONTSeq-Runtime, BAM-Speicher, Ausgabeverzeichnis und Reference-Lock.";
            await _launcher.VerifyPrerequisitesAsync(_settings, allowedRoot, referenceLock, cancellationToken);

            BackendStateText.Text = "Backend startet…";
            _launcher.Start(_settings, allowedRoot, referenceLock);
            _client = new OntSeqServiceClient(_settings.Port);
            var config = await _client.BootstrapAsync(
                TimeSpan.FromSeconds(30),
                () => _launcher.HasExited,
                () => _launcher.DiagnosticLog,
                cancellationToken);
            BackendStateText.Text = $"ONTSeq {config.Version} · lokal verbunden";

            _currentRunId = $"RUN_{DateTime.Now:yyyyMMdd_HHmmss}";
            _currentSampleId = sampleId;
            var request = new RunStartRequest(
                bam,
                sampleId,
                _currentRunId,
                genomeBuild,
                assay,
                assay == "adaptive_sampling" ? _settings.AdaptiveTargetBedWsl : null,
                assay == "adaptive_sampling" ? _settings.AdaptiveTargetBedVersion : null);

            if (assay == "adaptive_sampling" &&
                (string.IsNullOrWhiteSpace(request.TargetBed) || string.IsNullOrWhiteSpace(request.TargetBedVersion)))
            {
                throw new InvalidOperationException(
                    "Adaptive Sampling benötigt das freigegebene Analyse-ROI-BED und dessen Version. " +
                    "Diese Ressource ist noch nicht konfiguriert; der Lauf wird fail-closed gestoppt.");
            }

            await _client.StartRunAsync(request, cancellationToken);
            RunStateText.Text = "RUNNING";
            DetailText.Text = $"Analyse {_currentRunId} läuft. Die Bioinformatik arbeitet in WSL; dieses Fenster liest den geprüften Laufstatus.";
            await PollUntilFinishedAsync(_currentRunId, sampleId, cancellationToken);
        }
        catch (OperationCanceledException)
        {
            RunStateText.Text = "Unterbrochen";
            DetailText.Text = "Die Desktop-Überwachung wurde beendet. Das Backend wird beim Schließen der Anwendung gestoppt.";
        }
        catch (Exception error)
        {
            RunProgress.IsIndeterminate = false;
            RunStateText.Text = "FEHLER";
            DetailText.Text = error.Message;
            BackendStateText.Text = _launcher?.HasExited == true ? "Backend beendet" : BackendStateText.Text;
        }
        finally
        {
            StartButton.IsEnabled = true;
            BrowseButton.IsEnabled = true;
            SetupButton.IsEnabled = true;
        }
    }

    private async Task PollUntilFinishedAsync(string runId, string sampleId, CancellationToken cancellationToken)
    {
        if (_client is null) throw new InvalidOperationException("Backend-Client fehlt.");
        while (true)
        {
            cancellationToken.ThrowIfCancellationRequested();

            var persisted = await ProvenanceReader.ReadStagesAsync(
                _settings.OutputDirectoryWindows, runId, sampleId, cancellationToken);
            if (persisted.Count > 0) UpdateStages(persisted);

            var job = await _client.GetRunAsync(runId, cancellationToken);
            if (persisted.Count == 0 && job.Stages.Count > 0) UpdateStages(job.Stages);
            RunStateText.Text = job.State.ToUpperInvariant();

            if (!job.State.Equals("running", StringComparison.OrdinalIgnoreCase))
            {
                RunProgress.IsIndeterminate = false;
                if (job.Stages.Count > 0) UpdateStages(job.Stages);
                EnableOutputs(runId, sampleId);
                DetailText.Text = job.Detail;
                return;
            }

            await Task.Delay(1000, cancellationToken);
        }
    }

    private void UpdateStages(IReadOnlyCollection<StageSnapshot> stages)
    {
        StageItems.Clear();
        foreach (var stage in stages)
            StageItems.Add(new StageDisplay(stage.Title, stage.Status, stage.Reason));

        var concluded = stages.Count(s => s.Status is "COMPLETED" or "NO_CALL" or "FAILED" or "NOT_RUN");
        RunProgress.IsIndeterminate = false;
        RunProgress.Value = stages.Count == 0 ? 0 : 100.0 * concluded / stages.Count;
    }

    private bool TryValidateForm(out string bam, out string sampleId, out string build, out string assay)
    {
        bam = BamPathTextBox.Text.Trim();
        sampleId = SampleIdTextBox.Text.Trim();
        build = ((ComboBoxItem)GenomeBuildCombo.SelectedItem).Tag?.ToString() ?? "GRCh38";
        assay = ((ComboBoxItem)AssayCombo.SelectedItem).Tag?.ToString() ?? "lcwgs";

        if (!File.Exists(bam) || !bam.EndsWith(".bam", StringComparison.OrdinalIgnoreCase))
        {
            MessageBox.Show(this, "Bitte eine vorhandene BAM-Datei auswählen.", "ONTSeq", MessageBoxButton.OK, MessageBoxImage.Warning);
            return false;
        }
        if (BamIndexLocator.Find(bam) is null)
        {
            var shortIndex = Path.ChangeExtension(bam, ".bai");
            MessageBox.Show(this,
                "Der zugehörige BAM-Index fehlt. Erwartet wird entweder:\n" +
                bam + ".bai\noder:\n" + shortIndex,
                "ONTSeq", MessageBoxButton.OK, MessageBoxImage.Warning);
            return false;
        }
        if (!SafeId.IsMatch(sampleId))
        {
            MessageBox.Show(this,
                "Die Sample-ID muss 3–64 Zeichen lang sein und darf nur Buchstaben, Zahlen, Punkt, Unterstrich und Bindestrich enthalten.",
                "ONTSeq", MessageBoxButton.OK, MessageBoxImage.Warning);
            return false;
        }
        return true;
    }

    private void EnableOutputs(string runId, string sampleId)
    {
        var folder = ResultFolder(runId, sampleId);
        OpenFolderButton.IsEnabled = Directory.Exists(folder);
        OpenReportButton.IsEnabled = File.Exists(Path.Combine(folder, "reports", $"{sampleId}.report.html"));
        OpenExcelButton.IsEnabled = File.Exists(Path.Combine(folder, "reports", $"{sampleId}.results.xlsx"));
    }

    private string ResultFolder(string runId, string sampleId) =>
        Path.Combine(_settings.OutputDirectoryWindows, runId, sampleId);

    private void OpenReport_Click(object sender, RoutedEventArgs e) => OpenResultFile("report.html");
    private void OpenExcel_Click(object sender, RoutedEventArgs e) => OpenResultFile("results.xlsx");

    private void OpenResultFile(string suffix)
    {
        if (_currentRunId is null || _currentSampleId is null) return;
        var path = Path.Combine(ResultFolder(_currentRunId, _currentSampleId), "reports", $"{_currentSampleId}.{suffix}");
        if (File.Exists(path)) Process.Start(new ProcessStartInfo(path) { UseShellExecute = true });
    }

    private void OpenFolder_Click(object sender, RoutedEventArgs e)
    {
        if (_currentRunId is null || _currentSampleId is null) return;
        var folder = ResultFolder(_currentRunId, _currentSampleId);
        if (Directory.Exists(folder)) Process.Start(new ProcessStartInfo("explorer.exe", folder) { UseShellExecute = true });
    }

    private void SetPlaceholders()
    {
        StageItems.Clear();
        foreach (var title in new[] { "Inputprüfung", "Quality Control", "CNV", "Structural Variants", "Fusion", "Annotation", "Report" })
            StageItems.Add(new StageDisplay(title, "PENDING", ""));
    }

    private static string SanitizeSuggestedId(string value)
    {
        var cleaned = Regex.Replace(value, "[^A-Za-z0-9._-]", "_");
        if (cleaned.Length > 64) cleaned = cleaned[..64];
        if (cleaned.Length < 3) cleaned = "SAMPLE_" + cleaned;
        return cleaned;
    }

    private async void Window_Closed(object? sender, EventArgs e)
    {
        _analysisCts?.Cancel();
        _client?.Dispose();
        if (_launcher is not null) await _launcher.DisposeAsync();
        _analysisCts?.Dispose();
    }
}
