using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using System.Windows;
using Microsoft.Win32;

namespace ONTSeq.Desktop;

public partial class SetupWindow : Window
{
    private readonly DesktopSettings _settings;
    private readonly WslServiceLauncher _launcher = new();
    private CancellationTokenSource? _cts;

    public SetupWindow(DesktopSettings settings)
    {
        InitializeComponent();
        _settings = settings;
        _settings.ApplyProfileDefaults();
        ResourceRootTextBox.Text = _settings.ResourceRootWsl;
        SettingsPathText.Text = "Konfiguration: " + DesktopSettings.UserSettingsPath;
        RemoveAdaptiveBedButton.IsEnabled = _settings.HasAdaptiveTargetBedConfiguration;
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
        await RunBusyAsync(async token =>
        {
            var wsl = await _launcher.CheckWslAsync(_settings, token);
            WslStatusText.Text = Prefix(wsl.Ok) + wsl.Detail;

            if (!wsl.Ok)
            {
                BackendStatusText.Text = "— WSL muss zuerst funktionieren.";
                BundleStatusText.Text = "— WSL muss zuerst funktionieren.";
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
                await RefreshBundleStatusAsync(token);
            else
                BundleStatusText.Text = "✕ Runtime mit Resource-Registry erforderlich";
            await RefreshReferenceAsync("GRCh38", Grch38StatusText, token);
            await RefreshAdaptiveBedAsync(token);

            DetailText.Text = backend.Ok
                ? "System ist grundsätzlich bereit. Für neue Läufe muss der GRCh38-Bundle-Status vollständig sein; Legacy-Pfade werden nicht mit Profil-Bundles gemischt."
                : "Das Linux-Backend fehlt. Nutze 'Runtime installieren'; danach erneut prüfen.";
        });
    }

    private async Task RefreshBundleStatusAsync(CancellationToken token)
    {
        var status = await _launcher.CheckResourceBundlesAsync(_settings, token);
        BundleStatusText.Text = Prefix(status.Ok) +
            (string.IsNullOrWhiteSpace(status.Detail)
                ? $"Resource-Root geprüft: {_settings.ResourceRootWsl}"
                : status.Detail);
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
        await RunBusyAsync(async token =>
        {
            DetailText.Text = "Installiere gepinnte ONTSeq-Linux-Runtime mit QDNAseq/ACE in WSL…";
            var target = await _launcher.InstallBundledRuntimeAsync(_settings, archive, token);
            BackendStatusText.Text = "✓ Installiert: " + target;
            DetailText.Text = "ONTSeq Runtime wurde installiert und als Backend gespeichert.";
        });
        await RefreshAsync();
    }

    private async void ConfigureGrch38_Click(object sender, RoutedEventArgs e) =>
        await ConfigureReferenceAsync("GRCh38");

    private async void SaveResourceRoot_Click(object sender, RoutedEventArgs e)
    {
        await RunBusyAsync(token =>
        {
            token.ThrowIfCancellationRequested();
            _settings.ResourceRootWsl = DesktopSettings.NormalizeResourceRootWsl(
                ResourceRootTextBox.Text);
            _settings.SaveUserSettings();
            ResourceRootTextBox.Text = _settings.ResourceRootWsl;
            DetailText.Text = "Resource-Root gespeichert: " + _settings.ResourceRootWsl;
            return Task.CompletedTask;
        });
        await RefreshAsync();
    }

    private async void InstallBundle_Click(object sender, RoutedEventArgs e)
    {
        await RunBusyAsync(async token =>
        {
            _settings.ResourceRootWsl = DesktopSettings.NormalizeResourceRootWsl(
                ResourceRootTextBox.Text);
            _settings.SaveUserSettings();
            BundleStatusText.Text = "● GRCh38-Profilressourcen werden installiert…";
            DetailText.Text =
                "Installiere Referenz, HEMATOLOGY_v3, AML_AS_111_GRCh38_v1 und Profile nach " +
                $"{_settings.ResourceRootWsl}.";
            var detail = await _launcher.InstallGrch38ProfileResourcesAsync(_settings, token);
            BundleStatusText.Text = "✓ Installation abgeschlossen";
            DetailText.Text = string.IsNullOrWhiteSpace(detail)
                ? "Die vollständige GRCh38-Ressourcenfamilie wurde installiert."
                : detail;
        });
        await RefreshAsync();
    }

    private async void RepairBundle_Click(object sender, RoutedEventArgs e)
    {
        await RunBusyAsync(async token =>
        {
            _settings.ResourceRootWsl = DesktopSettings.NormalizeResourceRootWsl(
                ResourceRootTextBox.Text);
            _settings.SaveUserSettings();
            BundleStatusText.Text = "● Vollständige GRCh38-Ressourcenfamilie wird repariert…";
            DetailText.Text =
                "Repariere Referenz, HEMATOLOGY_v3, AML_AS_111_GRCh38_v1 und Profile unter " +
                $"{_settings.ResourceRootWsl}; manuelles Löschen ist nicht erforderlich.";
            var detail = await _launcher.RepairGrch38ProfileResourcesAsync(_settings, token);
            BundleStatusText.Text = "✓ Reparatur abgeschlossen";
            DetailText.Text = string.IsNullOrWhiteSpace(detail)
                ? "Die vollständige GRCh38-Ressourcenfamilie wurde repariert."
                : detail;
        });
        await RefreshAsync();
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

        await RunBusyAsync(async token =>
        {
            DetailText.Text = $"Erzeuge {build}-Reference-Lock…";
            var lockPath = await _launcher.ConfigureReferenceAsync(_settings, dialog.FileName, build, token);
            DetailText.Text = $"{build}-Reference-Lock erstellt: {lockPath}";
        });
        await RefreshAsync();
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

    private async Task RunBusyAsync(Func<CancellationToken, Task> action)
    {
        _cts?.Cancel();
        _cts?.Dispose();
        _cts = new CancellationTokenSource();
        SetBusy(true);
        try
        {
            await action(_cts.Token);
        }
        catch (OperationCanceledException)
        {
            DetailText.Text = "Vorgang abgebrochen.";
        }
        catch (Exception error)
        {
            DetailText.Text = error.Message;
            MessageBox.Show(this, error.Message, "ONTSeq Einrichtung", MessageBoxButton.OK, MessageBoxImage.Error);
        }
        finally
        {
            SetBusy(false);
        }
    }

    private void SetBusy(bool busy)
    {
        CheckButton.IsEnabled = !busy;
        InstallRuntimeButton.IsEnabled = !busy;
        SaveResourceRootButton.IsEnabled = !busy;
        InstallBundleButton.IsEnabled = !busy;
        RepairBundleButton.IsEnabled = !busy;
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
