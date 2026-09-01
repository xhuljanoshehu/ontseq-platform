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
    private bool _loadingSettings = true;

    public ObservableCollection<StageDisplay> StageItems { get; } = [];

    public MainWindow()
    {
        InitializeComponent();
        ProfileCombo.ItemsSource = DesktopProfiles.Supported;
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
            _loadingSettings = true;
            _settings = DesktopSettings.Load();
            SelectProfile(_settings.DefaultProfile);
            BackendStateText.Text =
                $"WSL: {_settings.WslDistribution} · Ressourcen: {_settings.ResourceRootWsl}";
            DetailText.Text =
                "ONTSeq löst FASTA, Annotation, Knowledge und gegebenenfalls Panel " +
                "automatisch aus dem gewählten GRCh38-Profil auf. Bundle-Status unter 'System einrichten' prüfen.";
            StartButton.IsEnabled = true;
        }
        catch (Exception error)
        {
            BackendStateText.Text = "Konfigurationsfehler";
            DetailText.Text = error.Message;
            StartButton.IsEnabled = false;
        }
        finally
        {
            _loadingSettings = false;
        }
    }

    private void SelectProfile(string profileId)
    {
        foreach (var item in ProfileCombo.Items.OfType<DesktopAnalysisProfile>())
        {
            if (!string.Equals(item.ProfileId, profileId, StringComparison.Ordinal)) continue;
            ProfileCombo.SelectedItem = item;
            return;
        }
        ProfileCombo.SelectedIndex = 0;
    }

    private void ProfileCombo_SelectionChanged(object sender, SelectionChangedEventArgs e)
    {
        if (ProfileCombo.SelectedItem is not DesktopAnalysisProfile profile) return;
        UpdateDetectedBuildText(profile);
        if (_loadingSettings) return;
        _settings.DefaultProfile = profile.ProfileId;
        try
        {
            _settings.SaveUserSettings();
        }
        catch (Exception error)
        {
            DetailText.Text = "Profil konnte nicht als Standard gespeichert werden: " + error.Message;
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
        if (ProfileCombo.SelectedItem is DesktopAnalysisProfile profile)
            UpdateDetectedBuildText(profile);
    }

    private void UpdateDetectedBuildText(DesktopAnalysisProfile profile)
    {
        DetectedBuildText.Text =
            $"{profile.GenomeBuild} · {profile.DictionaryLabel} · " +
            "BAM-Dictionary wird beim Start automatisch geprüft";
    }

    private async void Start_Click(object sender, RoutedEventArgs e)
    {
        if (!TryValidateForm(out var bam, out var sampleId, out var profile)) return;

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

            BackendStateText.Text = "WSL, Backend und GRCh38-Bundles werden geprüft…";
            DetectedBuildText.Text = "GRCh38 · BAM-Dictionary-Prüfung läuft im Backend…";
            DetailText.Text =
                "Vorprüfung: WSL2, ONTSeq-Runtime, BAM-Speicher, Resource-Root sowie Referenz-, Knowledge- und Panel-Bundles.";
            await _launcher.VerifyProfilePrerequisitesAsync(
                _settings,
                allowedRoot,
                profile.ProfileId,
                cancellationToken);

            BackendStateText.Text = "Backend startet…";
            _launcher.StartProfile(_settings, allowedRoot, profile.ProfileId);
            _client = new OntSeqServiceClient(_settings.Port);
            var config = await _client.BootstrapAsync(
                TimeSpan.FromSeconds(30),
                () => _launcher.HasExited,
                () => _launcher.DiagnosticLog,
                cancellationToken);
            if (config.Profiles is null || !config.Profiles.Contains(
                    profile.ProfileId, StringComparer.Ordinal))
            {
                throw new InvalidOperationException(
                    $"Das Backend kann das Profil {profile.ProfileId} aus dem Resource-Root " +
                    "nicht vollständig auflösen. Bitte in 'System einrichten' die vollständige " +
                    "Profilressourcen-Reparatur ausführen.");
            }
            BackendStateText.Text = $"ONTSeq {config.Version} · lokal verbunden";

            _currentSampleId = sampleId;
            var request = new RunStartRequest(
                bam,
                sampleId,
                null,
                profile.ProfileId,
                profile.GenomeBuild,
                profile.Assay);

            var started = await _client.StartRunAsync(request, cancellationToken);
            var runId = started.RunId;
            _currentRunId = runId;
            RunStateText.Text = "RUNNING";
            DetailText.Text =
                $"Analyse {runId} läuft mit {profile.ProfileId}. " +
                "Die Bioinformatik arbeitet in WSL; dieses Fenster liest den geprüften Laufstatus.";
            await PollUntilFinishedAsync(runId, sampleId, cancellationToken);
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
            var persistedBuild = await ProvenanceReader.ReadGenomeBuildAsync(
                _settings.OutputDirectoryWindows, runId, sampleId, cancellationToken);
            if (!string.IsNullOrWhiteSpace(persistedBuild))
                DetectedBuildText.Text = persistedBuild + " · automatisch geprüft und in Provenienz gespeichert";

            var job = await _client.GetRunAsync(runId, cancellationToken);
            if (!string.IsNullOrWhiteSpace(job.DetectedGenomeBuild))
                DetectedBuildText.Text = job.DetectedGenomeBuild + " · automatisch aus BAM-Dictionary erkannt";
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

    private bool TryValidateForm(
        out string bam,
        out string sampleId,
        out DesktopAnalysisProfile profile)
    {
        bam = BamPathTextBox.Text.Trim();
        sampleId = SampleIdTextBox.Text.Trim();
        var profileId = (ProfileCombo.SelectedItem as DesktopAnalysisProfile)?.ProfileId
                        ?? DesktopProfiles.DefaultProfileId;
        profile = DesktopProfiles.Require(profileId);

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
