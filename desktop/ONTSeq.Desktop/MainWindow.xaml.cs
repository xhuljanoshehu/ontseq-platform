using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using Microsoft.Win32;

namespace ONTSeq.Desktop;

public partial class MainWindow : Window
{
    private static readonly Regex SafeId = new("^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$", RegexOptions.Compiled);

    private DesktopSettings _settings = new();
    private WslServiceLauncher? _launcher;
    private OntSeqServiceClient? _client;
    private CancellationTokenSource? _analysisCts;
    private CancellationTokenSource? _workspaceCts;
    private CancellationTokenSource? _probeCts;
    private readonly MethylationProbeGuard _probeGuard = new();
    private readonly SemaphoreSlim _serviceGate = new(1, 1);
    private Task _probeTask = Task.CompletedTask;
    private int _probeOperations;
    private bool _desktopBusy;
    private bool _probeIsThorough;
    private bool _closing;
    private string? _probeCleanupWarning;
    private string? _serviceAllowedRootWindows;
    private string? _serviceSettingsKey;
    private string? _serviceInstanceId;
    private int? _servicePort;
    private string? _currentRunId;
    private string? _currentSampleId;
    private bool _loadingSettings = true;
    private IReadOnlyDictionary<string, ResourceFamilyState> _resourceFamilies =
        DesktopResourcePolicy.UnavailableFamilies("Ressourcenstatus wird geprüft.");

    public ObservableCollection<StageDisplay> StageItems { get; } = [];

    public MainWindow()
    {
        InitializeComponent();
        ProfileCombo.ItemsSource = DesktopResourcePolicy.BuildProfileOptions(_resourceFamilies);
        DataContext = this;
        SetPlaceholders();
    }

    private async void Window_Loaded(object sender, RoutedEventArgs e)
    {
        await ReloadSettingsStateAsync();
    }

    private async Task ReloadSettingsStateAsync()
    {
        InvalidateProbe("Konfiguration wird geprüft. Methylierung bitte anschließend erneut prüfen.");
        await _probeTask;
        var configuredProfile = DesktopProfiles.DefaultProfileId;
        try
        {
            _loadingSettings = true;
            _settings = DesktopSettings.Load();
            configuredProfile = _settings.DefaultProfile;
            _resourceFamilies = DesktopResourcePolicy.UnavailableFamilies(
                "Ressourcenstatus wird geprüft.");
            RebuildProfileOptions(configuredProfile);
            ApplySelectedProfileAvailability();
            BackendStateText.Text =
                $"WSL: {_settings.WslDistribution} · Ressourcen: {_settings.ResourceRootWsl}";

            _launcher ??= new WslServiceLauncher();
            _resourceFamilies = await _launcher.CheckResourceFamiliesAsync(
                _settings, CancellationToken.None);
            var decision = DesktopResourcePolicy.ResolveInitialProfile(
                configuredProfile, _resourceFamilies);
            var selectedProfileId = decision.Profile?.ProfileId ?? configuredProfile;
            RebuildProfileOptions(selectedProfileId);
            string? persistenceWarning = null;
            if (decision.IsFallback && decision.Profile is not null)
            {
                _settings.DefaultProfile = decision.Profile.ProfileId;
                try
                {
                    _settings.SaveUserSettings();
                }
                catch (Exception error)
                {
                    persistenceWarning =
                        "Die bereite Ersatzauswahl konnte nicht als Standard gespeichert werden: " +
                        error.Message;
                }
            }
            BackendStateText.Text = string.Join(" · ", DesktopResourcePolicy.GenomeBuilds.Select(
                build => $"{build}: {_resourceFamilies[build].AvailabilityLabel}"));
            DetailText.Text = string.Join(" ", new[]
            {
                decision.Notice,
                persistenceWarning,
                decision.Notice is null && persistenceWarning is null
                    ? "ONTSeq löst FASTA, Annotation, Knowledge und gegebenenfalls Panel " +
                      "automatisch aus dem gewählten Build-Profil auf."
                    : null
            }.Where(message => !string.IsNullOrWhiteSpace(message)));
        }
        catch (Exception error)
        {
            _resourceFamilies = DesktopResourcePolicy.UnavailableFamilies(error.Message);
            RebuildProfileOptions(configuredProfile);
            BackendStateText.Text = "Ressourcenstatus nicht verfügbar";
            DetailText.Text = error.Message;
        }
        finally
        {
            _loadingSettings = false;
            ApplySelectedProfileAvailability();
        }
    }

    private void RebuildProfileOptions(string profileId)
    {
        ProfileCombo.ItemsSource = DesktopResourcePolicy.BuildProfileOptions(_resourceFamilies);
        foreach (var item in ProfileCombo.Items.OfType<DesktopProfileOption>())
        {
            if (!string.Equals(item.Profile.ProfileId, profileId, StringComparison.Ordinal)) continue;
            ProfileCombo.SelectedItem = item;
            return;
        }
        ProfileCombo.SelectedIndex = 0;
    }

    private void ProfileCombo_SelectionChanged(object sender, SelectionChangedEventArgs e)
    {
        if (ProfileCombo.SelectedItem is not DesktopProfileOption option) return;
        UpdateDetectedBuildText(option.Profile);
        ApplySelectedProfileAvailability();
        if (_loadingSettings) return;
        InvalidateProbe("Profil geändert. Methylierung bitte erneut prüfen; vor dem Analysestart wird frisch geprüft.");
        if (!DesktopResourcePolicy.CanPersistProfile(option.Profile.ProfileId, _resourceFamilies))
            return;
        _settings.DefaultProfile = option.Profile.ProfileId;
        try
        {
            _settings.SaveUserSettings();
        }
        catch (Exception error)
        {
            DetailText.Text = "Profil konnte nicht als Standard gespeichert werden: " + error.Message;
        }
    }

    private void ApplySelectedProfileAvailability()
    {
        if (ProfileCombo.SelectedItem is not DesktopProfileOption option)
        {
            StartButton.IsEnabled = false;
            LiveWorkspaceButton.IsEnabled = false;
            ProfileStatusText.Text = "— Kein Analyseprofil ausgewählt.";
            ProfileStatusText.Foreground = Brushes.SlateGray;
            RefreshProbeActions();
            return;
        }

        var state = option.ResourceFamily;
        StartButton.IsEnabled = state.CanAnalyze && !_desktopBusy && _probeOperations == 0;
        LiveWorkspaceButton.IsEnabled = state.CanAnalyze && !_desktopBusy && _probeOperations == 0;
        RefreshProbeActions();
        var marker = state.Availability switch
        {
            ResourceFamilyAvailability.Ready => "✓",
            ResourceFamilyAvailability.NotInstalled => "○",
            ResourceFamilyAvailability.Incomplete => "!",
            _ => "—"
        };
        var action = state.Availability switch
        {
            ResourceFamilyAvailability.NotInstalled =>
                $" System einrichten → {state.GenomeBuild} → Installieren (Download).",
            ResourceFamilyAvailability.Incomplete =>
                $" System einrichten → {state.GenomeBuild} → Reparieren.",
            ResourceFamilyAvailability.Unavailable =>
                " System einrichten → System prüfen.",
            _ => string.Empty
        };
        ProfileStatusText.Text = $"{marker} {state.GenomeBuild}: {state.Detail}{action}";
        ProfileStatusText.Foreground = state.Availability switch
        {
            ResourceFamilyAvailability.Ready => Brushes.SeaGreen,
            ResourceFamilyAvailability.NotInstalled => Brushes.DarkGoldenrod,
            ResourceFamilyAvailability.Incomplete => Brushes.DarkOrange,
            _ => Brushes.SlateGray
        };
    }

    private async void Setup_Click(object sender, RoutedEventArgs e)
    {
        InvalidateProbe("Prüfung vor der Einrichtung beendet. Nach Änderungen bitte erneut prüfen.");
        await _probeTask;
        if (_launcher is { HasExited: false })
        {
            try
            {
                if (_client is null || (await _client.GetConfigAsync(CancellationToken.None)).Busy)
                    throw new InvalidOperationException(
                        "Das Backend arbeitet noch. Einrichtung erst nach Abschluss des laufenden Analyselaufs öffnen.");
            }
            catch (Exception error)
            {
                DetailText.Text = error.Message;
                return;
            }
        }
        var selectedBuild = (ProfileCombo.SelectedItem as DesktopProfileOption)?.Profile.GenomeBuild;
        var setup = new SetupWindow(_settings, selectedBuild) { Owner = this };
        setup.ShowDialog();
        await ReloadSettingsStateAsync();
    }

    private async void Browse_Click(object sender, RoutedEventArgs e)
    {
        var dialog = new OpenFileDialog
        {
            Title = "Oxford-Nanopore BAM auswählen",
            Filter = "BAM-Dateien (*.bam)|*.bam",
            CheckFileExists = true,
            Multiselect = false
        };
        if (dialog.ShowDialog(this) != true) return;
        InvalidateProbe("Neue BAM ausgewählt. Methylierungsinformationen werden geprüft…");
        BamPathTextBox.Text = dialog.FileName;
        var sample = Path.GetFileNameWithoutExtension(dialog.FileName);
        SampleIdTextBox.Text = SanitizeSuggestedId(sample);
        var index = BamIndexLocator.Find(dialog.FileName);
        DetailText.Text = index is not null
            ? $"BAM und BAM-Index gefunden ({Path.GetFileName(index)})."
            : "Hinweis: Erwartet wird <sample>.bam.bai oder <sample>.bai neben der BAM-Datei.";
        if (ProfileCombo.SelectedItem is DesktopProfileOption option)
            UpdateDetectedBuildText(option.Profile);
        ApplySelectedProfileAvailability();
        MethylationStateText.Text = "Methylierungsinformationen werden nach der Dateiauswahl geprüft…";
        if (index is null || ProfileCombo.SelectedItem is not DesktopProfileOption { IsEnabled: true } selected)
        {
            MethylationStateText.Text = "BAM-Index und installiertes Profil erforderlich. Methylierung wird vor dem Start geprüft.";
            return;
        }
        await BeginProbeAsync(dialog.FileName, selected.Profile, thorough: false);
    }

    private string SelectedProfileId =>
        (ProfileCombo.SelectedItem as DesktopProfileOption)?.Profile.ProfileId ?? "";

    private bool IsCurrentProbe(MethylationProbeScope scope) => !_closing &&
        _probeGuard.IsCurrent(scope, BamPathTextBox.Text.Trim(), SelectedProfileId, _serviceInstanceId);

    private void RefreshProbeActions()
    {
        if (ThoroughProbeButton is null) return;
        var ready = ProfileCombo.SelectedItem is DesktopProfileOption { IsEnabled: true };
        var probing = _probeOperations > 0;
        ThoroughProbeButton.IsEnabled = ready && !string.IsNullOrWhiteSpace(BamPathTextBox.Text) &&
            !_desktopBusy && !probing;
        CancelProbeButton.IsEnabled = !_desktopBusy && probing && _probeCts is { IsCancellationRequested: false };
        MethylationProgress.Visibility = probing ? Visibility.Visible : Visibility.Collapsed;
        SetupButton.IsEnabled = !_desktopBusy && !probing;
    }

    private void InvalidateProbe(string message)
    {
        _probeGuard.Invalidate();
        _probeCts?.Cancel();
        if (MethylationStateText is null) return;
        MethylationStateText.Text = message;
        MethylationExplanationText.Text = "Ein Prüfergebnis wählt keine Methylierungsauswertung automatisch aus.";
        MethylationProgressText.Text = "";
        MethylationDetailsText.Text = "";
        RefreshProbeActions();
    }

    private void ShowProbe(MethylationProbeResponse probe)
    {
        var presentation = MethylationProbePresentation.From(probe);
        MethylationStateText.Text = presentation.Title;
        MethylationExplanationText.Text = presentation.Explanation;
        MethylationProgressText.Text = presentation.Progress;
        MethylationDetailsText.Text = presentation.Detail;
    }

    private async void ThoroughProbe_Click(object sender, RoutedEventArgs e)
    {
        if (_desktopBusy || _probeOperations > 0) return;
        if (ProfileCombo.SelectedItem is not DesktopProfileOption { IsEnabled: true } option ||
            string.IsNullOrWhiteSpace(BamPathTextBox.Text)) return;
        await BeginProbeAsync(BamPathTextBox.Text.Trim(), option.Profile, thorough: true);
    }

    private Task BeginProbeAsync(string bam, DesktopAnalysisProfile profile, bool thorough)
    {
        var previous = _probeTask;
        InvalidateProbe(thorough ? "Gründliche Prüfung wird vorbereitet…" : "Methylierungs-Stichprobe wird geprüft…");
        _probeCleanupWarning = null;
        var cts = new CancellationTokenSource();
        _probeCts = cts;
        _probeIsThorough = thorough;
        var scope = _probeGuard.Capture(bam, profile.ProfileId, _serviceInstanceId);
        _probeTask = RunProbeAsync(previous, bam, profile, scope, thorough, cts);
        return _probeTask;
    }

    private async Task RunProbeAsync(Task previous, string bam, DesktopAnalysisProfile profile,
        MethylationProbeScope scope, bool thorough, CancellationTokenSource cts)
    {
        _probeOperations++;
        ApplySelectedProfileAvailability();
        _ = Dispatcher.BeginInvoke(new Action(() =>
        {
            if (!_closing && _probeGuard.IsCurrentSelection(scope, BamPathTextBox.Text.Trim(), SelectedProfileId))
                MethylationPanel.BringIntoView();
        }), System.Windows.Threading.DispatcherPriority.Loaded);
        OntSeqServiceClient? scanClient = null;
        string? scanId = null;
        var scanTerminal = false;
        var servicePrepared = false;
        var bamWsl = "";
        try
        {
            await previous; // Finish cleanup of the previous selection before opening another scan.
            cts.Token.ThrowIfCancellationRequested();
            bamWsl = PathBridge.WindowsToWsl(bam);
            await EnsureProfileServiceAsync(Path.GetDirectoryName(bam)!, profile, cts.Token);
            servicePrepared = true;
            if (!_probeGuard.IsCurrentSelection(scope, BamPathTextBox.Text.Trim(), SelectedProfileId)) return;
            scope = scope with { ServiceInstanceId = _serviceInstanceId };
            scanClient = _client!;
            if (!thorough)
            {
                var probe = await scanClient.ProbeMethylationAsync(bamWsl, cts.Token);
                if (IsCurrentProbe(scope)) ShowProbe(probe);
                return;
            }

            // Keep the short creation request alive even if selection changes. Its returned ID
            // is needed to cancel the server worker; every subsequent UI update is scope-checked.
            var snapshot = await scanClient.StartMethylationScanAsync(bamWsl, CancellationToken.None);
            scanId = snapshot.ScanId;
            while (true)
            {
                cts.Token.ThrowIfCancellationRequested();
                if (!IsCurrentProbe(scope)) return;
                MethylationProgressText.Text = MethylationProbePresentation.FormatProgress(
                    snapshot.CheckedReads, snapshot.ElapsedSeconds);
                if (snapshot.State != "running")
                {
                    scanTerminal = true;
                    ShowProbe(snapshot.Result ?? new MethylationProbeResponse(bamWsl, "unknown",
                        "Die Hintergrundprüfung lieferte kein abschließendes Prüfergebnis.",
                        snapshot.CheckedReads, false,
                        ReasonCode: snapshot.State == "cancelled" ? "cancelled" : "worker_failed",
                        ElapsedSeconds: snapshot.ElapsedSeconds, ScanMode: "thorough"));
                    return;
                }
                MethylationStateText.Text = "Gründliche Prüfung läuft";
                MethylationExplanationText.Text =
                    "Weitere Reads werden nach 5mC-Tags durchsucht. Bitte auf das Ergebnis warten oder die Prüfung abbrechen. Es wurde keine Analyse gestartet.";
                await Task.Delay(600, cts.Token);
                snapshot = await scanClient.GetMethylationScanAsync(scanId, bamWsl, cts.Token);
            }
        }
        catch (OperationCanceledException) when (cts.IsCancellationRequested) { }
        catch (Exception error)
        {
            if (IsCurrentProbe(scope)) ShowProbe(MethylationProbeFailures.FromException(
                bamWsl, error, thorough ? "thorough" : "quick", beforeProbe: !servicePrepared));
        }
        finally
        {
            if (scanClient is not null && scanId is not null && !scanTerminal)
            {
                try
                {
                    using var cleanupTimeout = new CancellationTokenSource(TimeSpan.FromSeconds(8));
                    var snapshot = await scanClient.CancelMethylationScanAsync(scanId, bamWsl, cleanupTimeout.Token);
                    while (snapshot.State == "running")
                    {
                        await Task.Delay(300, cleanupTimeout.Token);
                        snapshot = await scanClient.GetMethylationScanAsync(scanId, bamWsl, cleanupTimeout.Token);
                    }
                }
                catch (Exception error)
                {
                    _probeCleanupWarning = "Der Abbruch konnte noch nicht bestätigt werden. " + error.Message;
                }
            }
            if (ReferenceEquals(_probeCts, cts)) _probeCts = null;
            cts.Dispose();
            _probeOperations--;
            if (!_closing) ApplySelectedProfileAvailability();
        }
    }

    private async void CancelProbe_Click(object sender, RoutedEventArgs e)
    {
        if (_probeOperations == 0) return;
        var wasThorough = _probeIsThorough;
        var lastProgress = MethylationProgressText.Text;
        InvalidateProbe("Prüfung wird abgebrochen…");
        MethylationProgressText.Text = lastProgress;
        var scope = _probeGuard.Capture(BamPathTextBox.Text.Trim(), SelectedProfileId, _serviceInstanceId);
        var task = _probeTask;
        await task;
        if (!_probeGuard.IsCurrentSelection(scope, BamPathTextBox.Text.Trim(), SelectedProfileId) || _closing) return;
        MethylationStateText.Text = _probeCleanupWarning is not null ? "Abbruch noch nicht bestätigt" :
            wasThorough ? "Prüfung beendet; Ergebnis verworfen" : "Prüfanfrage abgebrochen";
        MethylationExplanationText.Text =
            "Es wurde keine Analyse gestartet. Der Methylierungsstatus bleibt offen; Sie können erneut prüfen oder vor dem Start nur die Genomanalyse wählen.";
        MethylationDetailsText.Text = _probeCleanupWarning ?? (wasThorough
            ? "Die Hintergrundprüfung ist beendet. Später eingehende Antworten dieser Prüfung werden verworfen."
            : "Die begrenzte Stichprobe kann im Dienst noch kurz auslaufen. Ihre Antwort wird nicht mehr übernommen.");
    }

    private void UpdateDetectedBuildText(DesktopAnalysisProfile profile)
    {
        DetectedBuildText.Text =
            $"{profile.GenomeBuild} · {profile.DictionaryLabel} · " +
            "BAM-Dictionary wird beim Start automatisch geprüft";
    }

    private async Task<ServiceConfigResponse> EnsureProfileServiceAsync(
        string allowedRoot, DesktopAnalysisProfile profile, CancellationToken cancellationToken)
    {
        await _serviceGate.WaitAsync(cancellationToken);
        try
        {
            return await EnsureProfileServiceCoreAsync(allowedRoot, profile, cancellationToken);
        }
        finally { _serviceGate.Release(); }
    }

    private async Task<ServiceConfigResponse> EnsureProfileServiceCoreAsync(
        string allowedRoot, DesktopAnalysisProfile profile, CancellationToken cancellationToken)
    {
        var settingsKey = ServiceSettingsKey();
        if (_launcher is { HasExited: false })
        {
            // Never terminate a service that may have been started/used through the browser.
            // Changing its filesystem scope requires an explicit Desktop restart.
            if (_client is null || _serviceSettingsKey != settingsKey ||
                !WslServiceLauncher.IsWithinAllowedRoot(allowedRoot, _serviceAllowedRootWindows))
                throw new InvalidOperationException(
                    "Der lokale Dienst läuft bereits mit einem anderen Eingabeordner oder einer anderen " +
                    "Konfiguration. Laufende Analysen zuerst abschließen, dann Desktop neu öffnen und " +
                    "die BAM vor dem Live-Workspace auswählen. Der Dienst wurde nicht beendet.");
            if (_serviceInstanceId is null)
                throw new InvalidOperationException(
                    "Die Identität des laufenden lokalen Dienstes ist nicht belegt. Desktop bitte neu öffnen.");
            return (await _client.GetConfigAsync(cancellationToken))
                .RequireLaunch(ServiceLaunchExpectationFor(
                    _serviceInstanceId, _serviceAllowedRootWindows!))
                .RequireProfile(profile.ProfileId);
        }

        DesktopOutputDirectory.EnsureExists(_settings.OutputDirectoryWindows);
        BackendStateText.Text = $"WSL, Backend und {profile.GenomeBuild}-Bundles werden geprüft…";
        await using (var verifier = new WslServiceLauncher())
        {
            await verifier.VerifyProfilePrerequisitesAsync(
                _settings, allowedRoot, profile.ProfileId, cancellationToken);
        }

        Exception? lastError = null;
        for (var attempt = 0; attempt < 2; attempt++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            var portSelection = attempt == 0
                ? WslServiceLauncher.SelectLoopbackPort(_settings.Port)
                : WslServiceLauncher.SelectEphemeralLoopbackPort();
            if (portSelection.UsedFallback)
            {
                BackendStateText.Text = attempt == 0
                    ? $"Port {_settings.Port} ist belegt · isolierte ONTSeq-Instanz startet auf " +
                      $"Port {portSelection.Port}"
                    : $"Portkollision sicher erkannt · neuer isolierter Start auf Port " +
                      $"{portSelection.Port}";
            }

            var instanceId = WslServiceLauncher.CreateServiceInstanceId();
            var expectation = ServiceLaunchExpectationFor(instanceId, allowedRoot);
            var candidateLauncher = new WslServiceLauncher();
            OntSeqServiceClient? candidateClient = null;
            try
            {
                candidateLauncher.StartProfile(
                    _settings, allowedRoot, profile.ProfileId, portSelection.Port, instanceId);
                candidateClient = new OntSeqServiceClient(portSelection.Port);
                var config = (await candidateClient.BootstrapAsync(
                        TimeSpan.FromSeconds(30), () => candidateLauncher.HasExited,
                        () => candidateLauncher.DiagnosticLog, expectation, cancellationToken))
                    .RequireProfile(profile.ProfileId);

                // Publish the candidate only after its launch identity and requested profile
                // have both been authenticated. Failed attempts never become reusable state.
                _client?.Dispose();
                if (_launcher is not null) await _launcher.DisposeAsync();
                _launcher = candidateLauncher;
                _client = candidateClient;
                _serviceAllowedRootWindows = allowedRoot;
                _serviceSettingsKey = settingsKey;
                _serviceInstanceId = instanceId;
                _servicePort = portSelection.Port;
                return config;
            }
            catch (Exception error)
            {
                lastError = error;
                candidateClient?.Dispose();
                await candidateLauncher.DisposeAsync();
                if (WslServiceLauncher.IsServiceLaunchCollision(
                        error, candidateLauncher.DiagnosticLog))
                {
                    if (attempt == 0) continue;
                    throw new InvalidOperationException(
                        "ONTSeq hat auch am neu gewählten Port keine eigene, eindeutig " +
                        "zuordenbare Dienstinstanz erhalten. Kein fremder Dienst wurde " +
                        "übernommen oder beendet.", error);
                }
                throw;
            }
        }

        throw new InvalidOperationException(
            "ONTSeq konnte nach einer erkannten lokalen Portkollision keine eigene " +
            "Dienstinstanz starten.", lastError);
    }

    private ServiceLaunchExpectation ServiceLaunchExpectationFor(
        string instanceId, string allowedRoot) => new(
            instanceId,
            _settings.ResourceRootWsl,
            PathBridge.WindowsToWsl(_settings.OutputDirectoryWindows),
            PathBridge.WindowsToWsl(allowedRoot));

    // The backend is intentionally multi-profile: the profile is bound to each POST /api/runs,
    // not to the service process. Its identity therefore contains launch scope/settings, while
    // RequireProfile above verifies the requested run capability on every reuse.
    private string ServiceSettingsKey() => string.Join("\n",
        _settings.WslDistribution, _settings.BackendCommand, _settings.RuntimeBinWsl,
        _settings.ModkitExecutableWsl,
        _settings.ResourceRootWsl, _settings.OutputDirectoryWindows, _settings.Port.ToString());

    private async void LiveWorkspace_Click(object sender, RoutedEventArgs e)
    {
        if (_desktopBusy || _probeOperations > 0) return;
        _desktopBusy = true;
        InvalidateProbe("Arbeitsplatz wird geöffnet. Methylierung wird vor einem Analysestart erneut geprüft.");
        LiveWorkspaceButton.IsEnabled = false;
        StartButton.IsEnabled = BrowseButton.IsEnabled = SetupButton.IsEnabled = false;
        ProfileCombo.IsEnabled = false;
        _workspaceCts?.Dispose();
        _workspaceCts = new CancellationTokenSource();
        try
        {
            var option = ProfileCombo.SelectedItem as DesktopProfileOption
                ?? throw new InvalidOperationException("Kein Analyseprofil ausgewählt.");
            if (!option.IsEnabled)
                throw new InvalidOperationException(
                    $"{option.Profile.GenomeBuild} ist im aktuellen Resource-Root nicht bereit. " +
                    "Bitte zuerst 'System einrichten' öffnen.");
            var allowedRoot = WslServiceLauncher.WorkspaceAllowedRootWindows(
                BamPathTextBox.Text.Trim(), _settings.OutputDirectoryWindows);
            var config = await EnsureProfileServiceAsync(
                allowedRoot, option.Profile, _workspaceCts.Token);
            BackendStateText.Text =
                $"ONTSeq {config.Version} · lokal verbunden · Port {_servicePort ?? _settings.Port}";
            DetailText.Text =
                $"Live-Workspace mit {option.Profile.ProfileId} geöffnet; kein Analyselauf gestartet. " +
                "Zugelassener Eingabeordner: " + allowedRoot +
                ". Desktop während der Nutzung geöffnet lassen.";
            Process.Start(new ProcessStartInfo(
                OntSeqServiceClient.WorkspaceUri(
                    _servicePort ?? _settings.Port, option.Profile.ProfileId).AbsoluteUri)
                { UseShellExecute = true });
        }
        catch (OperationCanceledException) { }
        catch (Exception error) { DetailText.Text = error.Message; }
        finally
        {
            _desktopBusy = false;
            BrowseButton.IsEnabled = SetupButton.IsEnabled = ProfileCombo.IsEnabled = true;
            ApplySelectedProfileAvailability();
        }
    }

    private async void Start_Click(object sender, RoutedEventArgs e)
    {
        if (_desktopBusy) return;
        if (_probeOperations > 0)
        {
            DetailText.Text = "Die Methylierungsprüfung läuft noch. Bitte das Ergebnis abwarten oder 'Prüfung abbrechen' wählen. Kein Analyselauf wurde gestartet.";
            return;
        }
        if (!TryValidateForm(out var bam, out var sampleId, out var profile)) return;

        _desktopBusy = true;
        InvalidateProbe("Methylierungsinformationen werden vor dem Start frisch geprüft…");

        StartButton.IsEnabled = false;
        BrowseButton.IsEnabled = false;
        SetupButton.IsEnabled = false;
        LiveWorkspaceButton.IsEnabled = false;
        ProfileCombo.IsEnabled = false;
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
        var servicePrepared = false;

        try
        {
            var allowedRoot = Path.GetDirectoryName(bam)
                              ?? throw new InvalidOperationException("BAM-Verzeichnis konnte nicht bestimmt werden.");

            DetectedBuildText.Text = $"{profile.GenomeBuild} · BAM-Dictionary-Prüfung läuft im Backend…";
            DetailText.Text =
                "Vorprüfung: WSL2, ONTSeq-Runtime, BAM-Speicher, Resource-Root sowie Referenz-, Knowledge- und Panel-Bundles.";
            var config = await EnsureProfileServiceAsync(allowedRoot, profile, cancellationToken);
            servicePrepared = true;
            BackendStateText.Text =
                $"ONTSeq {config.Version} · lokal verbunden · Port {_servicePort ?? _settings.Port}";
            if (config.Busy)
                throw new InvalidOperationException(
                    "Im lokalen Backend läuft bereits eine Analyse. Der bestehende Lauf wurde nicht verändert.");

            // Re-check the actual selected input immediately before submission. A prior
            // selection-time probe is informational and never grants opt-in for a run.
            MethylationStateText.Text = "Methylierungsinformationen werden vor dem Start geprüft…";
            MethylationProbeResponse probe;
            try
            {
                probe = await _client!.ProbeMethylationAsync(PathBridge.WindowsToWsl(bam), cancellationToken,
                    forceRefresh: true);
            }
            catch (Exception error) when (error is not OperationCanceledException || !cancellationToken.IsCancellationRequested)
            {
                probe = MethylationProbeFailures.FromException(PathBridge.WindowsToWsl(bam), error, "quick");
            }
            ShowProbe(probe);
            var probeExplanation = MethylationProbePresentation.From(probe).ConfirmationText;
            var includeMethylation = false;
            if (probe.Status == "detected" && probe.MethylationAvailable)
            {
                var decision = MessageBox.Show(this,
                    probeExplanation + "\n\n" +
                    "Soll zusätzlich zur Genomanalyse die regionale Methylierung ausgewertet werden?\n\n" +
                    "Ja: Methylierung ergänzen.\nNein: nur Genomanalyse.\nAbbrechen: noch keine Analyse starten.\n\n" +
                    "Die Methylierungsauswertung ist experimentell und nicht klinisch validiert.",
                    "Methylierung mitbeurteilen?", MessageBoxButton.YesNoCancel,
                    MessageBoxImage.Question, MessageBoxResult.Cancel);
                if (decision == MessageBoxResult.Cancel)
                {
                    RunStateText.Text = "Bereit";
                    DetailText.Text = "Kein Lauf gestartet. Der Analyseumfang kann erneut gewählt werden.";
                    return;
                }
                includeMethylation = decision == MessageBoxResult.Yes;
            }
            else if (probe.Status == "detected" && !probe.MethylationAvailable && MessageBox.Show(this,
                probeExplanation +
                "\n\nOhne Methylierung mit der Genomanalyse fortfahren?",
                "Methylierung noch nicht verfügbar", MessageBoxButton.OKCancel,
                MessageBoxImage.Information, MessageBoxResult.Cancel) != MessageBoxResult.OK)
            {
                RunStateText.Text = "Bereit";
                DetailText.Text = "Kein Lauf gestartet. Die Methylierungswerkzeuge müssen zunächst eingerichtet werden.";
                return;
            }
            else if (probe.Status == "unknown" && MessageBox.Show(this,
                probeExplanation + "\n\nOhne Methylierung mit der Genomanalyse fortfahren?\n" +
                "Mit Abbrechen können Sie zunächst 'Gründlicher prüfen' wählen.",
                "Methylierungsstatus unklar", MessageBoxButton.OKCancel,
                MessageBoxImage.Information, MessageBoxResult.Cancel) != MessageBoxResult.OK)
            {
                RunStateText.Text = "Bereit";
                DetailText.Text = "Kein Lauf gestartet. Bitte BAM-Prüfung oder Werkzeuge prüfen.";
                return;
            }
            MethylationStateText.Text = includeMethylation
                ? "Für diesen Lauf: Genomanalyse + regionale Methylierung."
                : "Für diesen Lauf: Genomanalyse ohne Methylierung.";
            _currentSampleId = sampleId;
            var request = new RunStartRequest(
                bam,
                sampleId,
                null,
                profile.ProfileId,
                profile.GenomeBuild,
                profile.Assay,
                IncludeMethylation: includeMethylation);

            var started = await _client!.StartRunAsync(request, cancellationToken);
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
            if (!servicePrepared)
                ShowProbe(MethylationProbeFailures.FromException(bam, error, "quick", beforeProbe: true));
            RunProgress.IsIndeterminate = false;
            RunStateText.Text = "FEHLER";
            DetailText.Text = error.Message;
            BackendStateText.Text = _launcher?.HasExited == true ? "Backend beendet" : BackendStateText.Text;
        }
        finally
        {
            _desktopBusy = false;
            RunProgress.IsIndeterminate = false;
            BrowseButton.IsEnabled = true;
            SetupButton.IsEnabled = true;
            ProfileCombo.IsEnabled = true;
            ApplySelectedProfileAvailability();
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
        var option = ProfileCombo.SelectedItem as DesktopProfileOption;
        profile = option?.Profile ?? DesktopProfiles.Require(DesktopProfiles.DefaultProfileId);

        if (option is null || !option.IsEnabled)
        {
            var build = option?.Profile.GenomeBuild ?? profile.GenomeBuild;
            MessageBox.Show(this,
                $"{build} ist im aktuellen Resource-Root nicht bereit. Bitte zuerst " +
                "'System einrichten' öffnen und die angezeigte Installations- oder Reparaturaktion ausführen.",
                "ONTSeq", MessageBoxButton.OK, MessageBoxImage.Warning);
            return false;
        }

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
        foreach (var title in new[] { "Inputprüfung", "Quality Control", "CNV", "Structural Variants", "Methylierung (optional)", "Fusion", "Annotation", "Report" })
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
        _closing = true;
        _analysisCts?.Cancel();
        _workspaceCts?.Cancel();
        _probeGuard.Invalidate();
        _probeCts?.Cancel();
        await _probeTask;
        _client?.Dispose();
        if (_launcher is not null) await _launcher.DisposeAsync();
        _analysisCts?.Dispose();
        _workspaceCts?.Dispose();
        _probeCts?.Dispose();
    }
}
