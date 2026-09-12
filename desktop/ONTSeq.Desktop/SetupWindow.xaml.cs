using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using System.Windows;
using System.Windows.Media;
using Microsoft.Win32;

namespace ONTSeq.Desktop;

public partial class SetupWindow : Window
{
    private readonly DesktopSettings _settings;
    private readonly WslServiceLauncher _launcher = new();
    private CancellationTokenSource? _cts;
    private IReadOnlyDictionary<string, ResourceFamilyState> _resourceFamilies =
        DesktopResourcePolicy.UnavailableFamilies("Noch nicht geprüft.");
    private bool _busy;

    public SetupWindow(DesktopSettings settings, string? initialResourceBuild = null)
    {
        InitializeComponent();
        _settings = settings;
        _settings.ApplyProfileDefaults();
        ResourceBuildCombo.ItemsSource = DesktopResourcePolicy.GenomeBuilds;
        var defaultBuild = DesktopProfiles.Require(_settings.DefaultProfile).GenomeBuild;
        ResourceBuildCombo.SelectedItem = DesktopResourcePolicy.GenomeBuilds.Contains(
            initialResourceBuild, StringComparer.Ordinal)
            ? initialResourceBuild
            : defaultBuild;
        ResourceRootTextBox.Text = _settings.ResourceRootWsl;
        SettingsPathText.Text = "Konfiguration: " + DesktopSettings.UserSettingsPath;
        RemoveAdaptiveBedButton.IsEnabled = _settings.HasAdaptiveTargetBedConfiguration;
        UpdateResourceStatusDisplay();
    }

    private async void Window_Loaded(object sender, RoutedEventArgs e)
    {
        await RefreshAsync();
    }

    private async void Check_Click(object sender, RoutedEventArgs e)
    {
        await RefreshAsync();
    }

    private async Task RefreshAsync()
    {
        SetResourceFamiliesUnavailable("Ressourcenstatus wird erneut geprüft.");
        await RunBusyAsync(async token =>
        {
            var wsl = await _launcher.CheckWslAsync(_settings, token);
            WslStatusText.Text = Prefix(wsl.Ok) + wsl.Detail;

            if (!wsl.Ok)
            {
                BackendStatusText.Text = "— WSL muss zuerst funktionieren.";
                SetResourceFamiliesUnavailable("WSL muss zuerst funktionieren.");
                Grch38StatusText.Text = _settings.TryReferenceLockFor("GRCh38", out _)
                    ? "— Legacy-Pfad gespeichert; WSL muss für die Prüfung funktionieren."
                    : "— Kein Legacy-Pfad gespeichert";
                AdaptiveBedStatusText.Text = _settings.HasAdaptiveTargetBedConfiguration
                    ? "— Legacy-Pfad gespeichert; WSL muss für die Prüfung funktionieren."
                    : "— Kein Legacy-Pfad gespeichert";
                SelfTestStatusText.Text = "— Nicht möglich, solange WSL fehlt.";
                return;
            }

            var backend = await _launcher.CheckBackendAsync(_settings, token);
            BackendStatusText.Text = Prefix(backend.Ok) + backend.Detail;
            if (backend.Ok)
                await RefreshBundleStatusesAsync(token);
            else
                SetResourceFamiliesUnavailable("Runtime mit Resource-Registry erforderlich.");
            await RefreshReferenceAsync("GRCh38", Grch38StatusText, token);
            await RefreshAdaptiveBedAsync(token);

            DetailText.Text = backend.Ok
                ? "System ist grundsätzlich bereit. Für neue Läufe muss die gewählte Build-Familie vollständig sein; die andere Familie ist optional. Legacy-Pfade werden nicht mit Profil-Bundles gemischt."
                : "Das Linux-Backend fehlt. Nutze 'Runtime installieren'; danach erneut prüfen.";
        });
    }

    private async Task RefreshBundleStatusesAsync(CancellationToken token)
    {
        _resourceFamilies = await _launcher.CheckResourceFamiliesAsync(_settings, token);
        UpdateResourceStatusDisplay();
    }

    private async Task RefreshAfterActionAsync(bool completed)
    {
        var failureDetail = completed ? null : DetailText.Text;
        await RefreshAsync();
        // A status refresh may discover partial installation, but must not replace the
        // operation's diagnostic with the generic "System ist grundsätzlich bereit".
        if (failureDetail is not null) DetailText.Text = failureDetail;
    }

    private string SelectedResourceBuild => ResourceBuildCombo.SelectedItem as string ?? "GRCh38";

    private void ResourceBuildCombo_SelectionChanged(
        object sender, System.Windows.Controls.SelectionChangedEventArgs e)
    {
        if (!IsLoaded) return;
        UpdateResourceActionButtons();
    }

    private void ResourceRootTextBox_TextChanged(
        object sender, System.Windows.Controls.TextChangedEventArgs e)
    {
        if (!IsLoaded || _busy) return;
        SetResourceFamiliesUnavailable(
            "Resource-Root geändert; mit 'Übernehmen' speichern und beide Build-Familien neu prüfen.");
    }

    private void SetResourceFamiliesUnavailable(string detail)
    {
        _resourceFamilies = DesktopResourcePolicy.UnavailableFamilies(detail);
        UpdateResourceStatusDisplay();
    }

    private void UpdateResourceStatusDisplay()
    {
        UpdateResourceStatusLine("GRCh37", Grch37BundleStatusText);
        UpdateResourceStatusLine("GRCh38", Grch38BundleStatusText);
        var readyBuilds = DesktopResourcePolicy.GenomeBuilds
            .Where(build => _resourceFamilies[build].CanAnalyze).ToArray();
        ResourceAvailabilityHintText.Text = readyBuilds.Length switch
        {
            2 => "Beide Referenzfamilien sind lokal bereit. Jeder Lauf verwendet ausschließlich " +
                 "die Familie des gewählten Analyseprofils.",
            1 => $"{readyBuilds[0]} kann bereits verwendet werden. Die andere Familie ist für " +
                 "diese Profile nicht erforderlich und kann zusätzlich eingerichtet werden.",
            _ => "Noch keine Referenzfamilie als bereit bestätigt. Für einen Lauf wird nur die " +
                 "Familie des gewählten Analyseprofils benötigt; beide können parallel bereitstehen."
        };
        UpdateResourceActionButtons();
    }

    private void UpdateResourceStatusLine(
        string genomeBuild,
        System.Windows.Controls.TextBlock target)
    {
        var state = _resourceFamilies[genomeBuild];
        var marker = state.Availability switch
        {
            ResourceFamilyAvailability.Ready => "✓",
            ResourceFamilyAvailability.NotInstalled => "○",
            ResourceFamilyAvailability.Incomplete => "!",
            _ => "—"
        };
        target.Text = $"{marker} {genomeBuild}: {state.Detail}";
        target.Foreground = state.Availability switch
        {
            ResourceFamilyAvailability.Ready => Brushes.SeaGreen,
            ResourceFamilyAvailability.NotInstalled => Brushes.DarkGoldenrod,
            ResourceFamilyAvailability.Incomplete => Brushes.DarkOrange,
            _ => Brushes.SlateGray
        };
    }

    private void SetResourceProgress(string genomeBuild, string text)
    {
        var target = genomeBuild == "GRCh37"
            ? Grch37BundleStatusText
            : Grch38BundleStatusText;
        target.Text = "● " + genomeBuild + ": " + text;
        target.Foreground = Brushes.SteelBlue;
    }

    private void UpdateResourceActionButtons()
    {
        var state = _resourceFamilies[SelectedResourceBuild];
        InstallBundleButton.IsEnabled = !_busy && state.CanInstall;
        RepairBundleButton.IsEnabled = !_busy && state.CanRepair;
        InstallBundleButton.Content = $"{state.GenomeBuild} installieren";
        RepairBundleButton.Content = state.Availability == ResourceFamilyAvailability.Ready
            ? $"{state.GenomeBuild} prüfen / reparieren"
            : $"{state.GenomeBuild} reparieren";
        InstallBundleButton.ToolTip = state.CanInstall
            ? $"{state.GenomeBuild} in diesem Resource-Root installieren; Internetzugang und mehrere GB Speicherplatz sind erforderlich."
            : "Installieren ist nur für eine noch nicht vorhandene Build-Familie verfügbar.";
        RepairBundleButton.ToolTip = state.CanRepair
            ? $"Vorhandene {state.GenomeBuild}-Ressourcen vollständig prüfen und fehlende oder beschädigte Dateien reparieren; ein Download kann erforderlich sein."
            : "Reparieren ist erst möglich, wenn die Referenzfamilie installiert ist.";
    }

    private async Task RefreshReferenceAsync(
        string build,
        System.Windows.Controls.TextBlock target,
        CancellationToken token)
    {
        if (!_settings.TryReferenceLockFor(build, out var referenceLock))
        {
            target.Text = "— Nicht konfiguriert";
            return;
        }

        var check = await _launcher.CheckReferenceAsync(
            _settings, referenceLock, build, token);
        target.Text = Prefix(check.Ok) + (string.IsNullOrWhiteSpace(check.Detail)
            ? $"Reference-Lock konnte nicht geprüft werden: {referenceLock}"
            : check.Detail);
    }

    private async Task RefreshAdaptiveBedAsync(CancellationToken token)
    {
        if (!_settings.HasAdaptiveTargetBedConfiguration)
        {
            AdaptiveBedStatusText.Text = "— Nicht konfiguriert";
            return;
        }

        if (string.IsNullOrWhiteSpace(_settings.AdaptiveTargetBedWsl) ||
            string.IsNullOrWhiteSpace(_settings.AdaptiveTargetBedVersion))
        {
            AdaptiveBedStatusText.Text =
                "✕ Unvollständige Konfiguration – bitte entfernen und das Analyse-BED neu wählen.";
            return;
        }

        var check = await WslServiceLauncher.RunWslAsync(
            _settings.WslDistribution,
            ["sh", "-lc", $"test -s {ShellQuote(_settings.AdaptiveTargetBedWsl)}"],
            token);
        AdaptiveBedStatusText.Text = check.ExitCode == 0
            ? $"✓ {_settings.AdaptiveTargetBedVersion}"
            : $"✕ Konfiguriert, aber in WSL nicht auffindbar: {_settings.AdaptiveTargetBedWsl}";
    }

    private async void InstallRuntime_Click(object sender, RoutedEventArgs e)
    {
        var archive = Path.Combine(AppContext.BaseDirectory, "runtime", "ontseq-linux-runtime.tar.gz");
        var completed = await RunBusyAsync(async token =>
        {
            DetailText.Text = "Prüfe Runtime-/Core-Prüfsummen und installiere in einen neuen WSL-Prefix; bestehende Runtime bleibt unverändert…";
            var target = await _launcher.InstallBundledRuntimeAsync(_settings, archive, token);
            BackendStatusText.Text = "✓ Installiert: " + target;
            DetailText.Text = "ONTSeq Runtime wurde installiert und als Backend gespeichert.";
        });
        await RefreshAfterActionAsync(completed);
    }

    private async void ConfigureGrch38_Click(object sender, RoutedEventArgs e) =>
        await ConfigureReferenceAsync("GRCh38");

    private async void SaveResourceRoot_Click(object sender, RoutedEventArgs e)
    {
        var completed = await RunBusyAsync(token =>
        {
            token.ThrowIfCancellationRequested();
            _settings.ResourceRootWsl = DesktopSettings.NormalizeResourceRootWsl(
                ResourceRootTextBox.Text);
            _settings.SaveUserSettings();
            ResourceRootTextBox.Text = _settings.ResourceRootWsl;
            DetailText.Text = "Resource-Root gespeichert: " + _settings.ResourceRootWsl;
            return Task.CompletedTask;
        });
        await RefreshAfterActionAsync(completed);
    }

    private async void InstallBundle_Click(object sender, RoutedEventArgs e)
    {
        var build = SelectedResourceBuild;
        var completed = await RunBusyAsync(async token =>
        {
            _settings.ResourceRootWsl = DesktopSettings.NormalizeResourceRootWsl(
                ResourceRootTextBox.Text);
            _settings.SaveUserSettings();
            SetResourceProgress(build, "Profilressourcen werden installiert…");
            DetailText.Text =
                "Installiere " + string.Join(", ", WslServiceLauncher.ManagedResourceBundleIds(build)) + " und Profile nach " +
                $"{_settings.ResourceRootWsl}.";
            var detail = await _launcher.InstallProfileResourcesAsync(_settings, build, token);
            SetResourceProgress(build, "Installation abgeschlossen; Status wird neu geprüft…");
            DetailText.Text = string.IsNullOrWhiteSpace(detail)
                ? $"Die vollständige {build}-Ressourcenfamilie wurde installiert."
                : detail;
        }, $"{build}: Profilressourcen konnten nicht installiert werden.");
        await RefreshAfterActionAsync(completed);
    }

    private async void RepairBundle_Click(object sender, RoutedEventArgs e)
    {
        var build = SelectedResourceBuild;
        var completed = await RunBusyAsync(async token =>
        {
            _settings.ResourceRootWsl = DesktopSettings.NormalizeResourceRootWsl(
                ResourceRootTextBox.Text);
            _settings.SaveUserSettings();
            SetResourceProgress(build, "Vollständige Ressourcenfamilie wird geprüft und repariert…");
            DetailText.Text =
                "Repariere " + string.Join(", ", WslServiceLauncher.ManagedResourceBundleIds(build)) + " und Profile unter " +
                $"{_settings.ResourceRootWsl}; manuelles Löschen ist nicht erforderlich.";
            var detail = await _launcher.RepairProfileResourcesAsync(_settings, build, token);
            SetResourceProgress(build, "Reparatur abgeschlossen; Status wird neu geprüft…");
            DetailText.Text = string.IsNullOrWhiteSpace(detail)
                ? $"Die vollständige {build}-Ressourcenfamilie wurde repariert."
                : detail;
        }, $"{build}: Profilressourcen konnten nicht repariert werden.");
        await RefreshAfterActionAsync(completed);
    }

    private async Task ConfigureReferenceAsync(string build)
    {
        var dialog = new OpenFileDialog
        {
            Title = $"{build}-Referenz auswählen",
            Filter = "Referenz FASTA/FAI (*.fa;*.fasta;*.fna;*.fai)|*.fa;*.fasta;*.fna;*.fai|Alle Dateien (*.*)|*.*",
            CheckFileExists = true,
            Multiselect = false
        };
        if (dialog.ShowDialog(this) != true) return;

        var confirm = MessageBox.Show(
            this,
            $"Die gewählte Datei muss exakt zu der Referenz passen, gegen die dein BAM ausgerichtet wurde.\n\n{dialog.FileName}\n\nAls {build} konfigurieren?",
            "ONTSeq Reference-Lock",
            MessageBoxButton.YesNo,
            MessageBoxImage.Warning);
        if (confirm != MessageBoxResult.Yes) return;

        var completed = await RunBusyAsync(async token =>
        {
            DetailText.Text = $"Erzeuge {build}-Reference-Lock…";
            var lockPath = await _launcher.ConfigureReferenceAsync(_settings, dialog.FileName, build, token);
            DetailText.Text = $"{build}-Reference-Lock erstellt: {lockPath}";
        });
        await RefreshAfterActionAsync(completed);
    }

    private async void ConfigureAdaptiveBed_Click(object sender, RoutedEventArgs e)
    {
        var dialog = new OpenFileDialog
        {
            Title = "Adaptive-Sampling Analyse-ROI-BED auswählen",
            Filter = "BED-Dateien (*.bed)|*.bed|Alle Dateien (*.*)|*.*",
            CheckFileExists = true,
            Multiselect = false
        };
        if (dialog.ShowDialog(this) != true) return;

        await RunBusyAsync(async token =>
        {
            ValidateBed(dialog.FileName);
            token.ThrowIfCancellationRequested();

            string hash;
            await using (var stream = new FileStream(
                             dialog.FileName,
                             FileMode.Open,
                             FileAccess.Read,
                             FileShare.Read,
                             1024 * 1024,
                             useAsync: true))
            {
                hash = Convert.ToHexString(await SHA256.HashDataAsync(stream, token)).ToLowerInvariant();
            }

            var resourceDirectory = Path.Combine(
                Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
                "ONTSeq", "resources", "adaptive_sampling");
            Directory.CreateDirectory(resourceDirectory);
            var destination = Path.Combine(resourceDirectory, $"analysis_roi.{hash[..12]}.bed");
            var temporary = destination + ".tmp";
            try
            {
                File.Copy(dialog.FileName, temporary, overwrite: true);
                File.Move(temporary, destination, overwrite: true);
            }
            finally
            {
                if (File.Exists(temporary)) File.Delete(temporary);
            }

            _settings.AdaptiveTargetBedWsl = PathBridge.WindowsToWsl(destination);
            _settings.AdaptiveTargetBedVersion =
                $"{Path.GetFileName(dialog.FileName)}@sha256:{hash}";
            _settings.SaveUserSettings();

            var check = await WslServiceLauncher.RunWslAsync(
                _settings.WslDistribution,
                ["sh", "-lc", $"test -s {ShellQuote(_settings.AdaptiveTargetBedWsl)}"],
                token);
            if (check.ExitCode != 0)
                throw new InvalidOperationException(
                    "Das gespeicherte Analyse-BED ist in WSL nicht erreichbar.\n" + check.StdErr);

            AdaptiveBedStatusText.Text = "✓ " + _settings.AdaptiveTargetBedVersion;
            DetailText.Text =
                "Adaptive-Sampling Analyse-ROI gespeichert und per SHA256 versioniert: " + destination;
        });
    }

    private void RemoveAdaptiveBed_Click(object sender, RoutedEventArgs e)
    {
        if (!_settings.HasAdaptiveTargetBedConfiguration) return;

        var configuredVersion = string.IsNullOrWhiteSpace(_settings.AdaptiveTargetBedVersion)
            ? "Unbekannte/teilweise Konfiguration"
            : _settings.AdaptiveTargetBedVersion;
        var confirm = MessageBox.Show(
            this,
            $"Adaptive-Sampling Analyse-BED aus der ONTSeq-Konfiguration entfernen?\n\n{configuredVersion}\n\n" +
            "Die ursprüngliche BED-Datei auf deinem Rechner wird nicht gelöscht.",
            "Adaptive Sampling ROI entfernen",
            MessageBoxButton.YesNo,
            MessageBoxImage.Question);
        if (confirm != MessageBoxResult.Yes) return;

        _settings.ClearAdaptiveTargetBed();
        _settings.SaveUserSettings();
        AdaptiveBedStatusText.Text = "— Nicht konfiguriert";
        RemoveAdaptiveBedButton.IsEnabled = false;
        DetailText.Text =
            "Adaptive-Sampling Analyse-ROI aus der ONTSeq-Konfiguration entfernt. " +
            "Die ursprüngliche BED-Datei wurde nicht gelöscht.";
    }

    private static void ValidateBed(string path)
    {
        var intervals = 0;
        var lineNumber = 0;
        foreach (var raw in File.ReadLines(path))
        {
            lineNumber++;
            var line = raw.Trim();
            if (line.Length == 0 || line.StartsWith('#')) continue;
            var fields = line.Split('\t');
            if (fields.Length < 3 || string.IsNullOrWhiteSpace(fields[0]))
                throw new InvalidDataException($"BED-Zeile {lineNumber} hat weniger als drei gültige Spalten.");
            if (!long.TryParse(fields[1], NumberStyles.None, CultureInfo.InvariantCulture, out var start) ||
                !long.TryParse(fields[2], NumberStyles.None, CultureInfo.InvariantCulture, out var end) ||
                start < 0 || end <= start)
            {
                throw new InvalidDataException(
                    $"BED-Zeile {lineNumber} hat ungültige 0-basierte Start/End-Koordinaten.");
            }
            intervals++;
        }
        if (intervals == 0)
            throw new InvalidDataException("Das gewählte BED enthält keine auswertbaren Intervalle.");
    }

    private async void SelfTest_Click(object sender, RoutedEventArgs e)
    {
        await RunBusyAsync(async token =>
        {
            SelfTestStatusText.Text = "● Läuft…";
            DetailText.Text =
                "Systemtest läuft: samtools, Cramino, Sniffles2, QDNAseq/ACE, kanonische Pipeline, HTML/Excel, Release-Checksummen und Resume werden mit deterministischen synthetischen Daten geprüft…";
            var folder = await _launcher.RunSelfTestAsync(_settings, token);
            SelfTestStatusText.Text = "✓ PASS";
            DetailText.Text =
                "Vollständiger Engineering-Systemtest erfolgreich. Ergebnis und Prüfnachweise: " + folder;
            if (MessageBox.Show(this, "Selbsttest erfolgreich. Ergebnisordner öffnen?", "ONTSeq", MessageBoxButton.YesNo, MessageBoxImage.Information)
                == MessageBoxResult.Yes)
            {
                Process.Start(new ProcessStartInfo("explorer.exe", folder) { UseShellExecute = true });
            }
        });
    }

    private async Task<bool> RunBusyAsync(
        Func<CancellationToken, Task> action, string? failureContext = null)
    {
        _cts?.Cancel();
        _cts?.Dispose();
        _cts = new CancellationTokenSource();
        SetBusy(true);
        try
        {
            await action(_cts.Token);
            return true;
        }
        catch (OperationCanceledException)
        {
            DetailText.Text = failureContext is null
                ? "Vorgang abgebrochen."
                : failureContext + "\nVorgang abgebrochen.";
            return false;
        }
        catch (Exception error)
        {
            DetailText.Text = failureContext is null
                ? error.Message
                : failureContext + "\n" + error.Message;
            MessageBox.Show(this, DetailText.Text, "ONTSeq Einrichtung", MessageBoxButton.OK, MessageBoxImage.Error);
            return false;
        }
        finally
        {
            SetBusy(false);
        }
    }

    private void SetBusy(bool busy)
    {
        _busy = busy;
        CheckButton.IsEnabled = !busy;
        InstallRuntimeButton.IsEnabled = !busy;
        SaveResourceRootButton.IsEnabled = !busy;
        ResourceBuildCombo.IsEnabled = !busy;
        UpdateResourceActionButtons();
        ConfigureGrch38Button.IsEnabled = !busy;
        ConfigureAdaptiveBedButton.IsEnabled = !busy;
        RemoveAdaptiveBedButton.IsEnabled = !busy && _settings.HasAdaptiveTargetBedConfiguration;
        SelfTestButton.IsEnabled = !busy;
    }

    private void Done_Click(object sender, RoutedEventArgs e)
    {
        _settings.SaveUserSettings();
        DialogResult = true;
        Close();
    }

    protected override void OnClosed(EventArgs e)
    {
        _cts?.Cancel();
        _cts?.Dispose();
        _launcher.DisposeAsync().AsTask().GetAwaiter().GetResult();
        base.OnClosed(e);
    }

    private static string Prefix(bool ok) => ok ? "✓ " : "✕ ";
    private static string ShellQuote(string value) => "'" + value.Replace("'", "'\"'\"'") + "'";
}
