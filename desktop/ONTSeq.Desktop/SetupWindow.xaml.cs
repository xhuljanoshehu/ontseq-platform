using System.Diagnostics;
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
        SettingsPathText.Text = "Konfiguration: " + DesktopSettings.UserSettingsPath;
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
                SelfTestStatusText.Text = "— Nicht möglich, solange WSL fehlt.";
                return;
            }

            var backend = await _launcher.CheckBackendAsync(_settings, token);
            BackendStatusText.Text = Prefix(backend.Ok) + backend.Detail;
            await RefreshReferenceAsync("GRCh38", Grch38StatusText, token);
            await RefreshReferenceAsync("GRCh37", Grch37StatusText, token);

            DetailText.Text = backend.Ok
                ? "System ist grundsätzlich bereit. Konfiguriere mindestens den Reference-Lock des BAM und führe danach den Selbsttest aus."
                : "Das Linux-Backend fehlt. Nutze 'Runtime installieren'; danach erneut prüfen.";
        });
    }

    private async Task RefreshReferenceAsync(string build, System.Windows.Controls.TextBlock target, CancellationToken token)
    {
        if (!_settings.TryReferenceLockFor(build, out var referenceLock))
        {
            target.Text = "— Nicht konfiguriert";
            return;
        }

        var check = await WslServiceLauncher.RunWslAsync(
            _settings.WslDistribution, ["sh", "-lc", $"test -f {ShellQuote(referenceLock)}"], token);
        target.Text = check.ExitCode == 0
            ? $"✓ Konfiguriert: {referenceLock}"
            : $"✕ Konfiguriert, aber in WSL nicht auffindbar: {referenceLock}";
    }

    private async void InstallRuntime_Click(object sender, RoutedEventArgs e)
    {
        var archive = Path.Combine(AppContext.BaseDirectory, "runtime", "ontseq-linux-runtime.tar.gz");
        await RunBusyAsync(async token =>
        {
            DetailText.Text = "Installiere gepinnte ONTSeq-Linux-Runtime in WSL…";
            var target = await _launcher.InstallBundledRuntimeAsync(_settings, archive, token);
            BackendStatusText.Text = "✓ Installiert: " + target;
            DetailText.Text = "ONTSeq Runtime wurde installiert und als Backend gespeichert.";
        });
        await RefreshAsync();
    }

    private async void ConfigureGrch38_Click(object sender, RoutedEventArgs e) =>
        await ConfigureReferenceAsync("GRCh38");

    private async void ConfigureGrch37_Click(object sender, RoutedEventArgs e) =>
        await ConfigureReferenceAsync("GRCh37");

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

    private async void SelfTest_Click(object sender, RoutedEventArgs e)
    {
        await RunBusyAsync(async token =>
        {
            SelfTestStatusText.Text = "● Läuft…";
            DetailText.Text = "Erzeuge synthetischen BAM und führe samtools, Cramino, Sniffles2 und Reporting aus…";
            var folder = await _launcher.RunSelfTestAsync(_settings, token);
            SelfTestStatusText.Text = "✓ PASS";
            DetailText.Text = "Engineering-Selbsttest erfolgreich. Ergebnis: " + folder;
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
