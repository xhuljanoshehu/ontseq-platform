param(
    [string]$OutputDir = "results\quick-test",
    [switch]$NoInstall,
    [switch]$SkipChecks,
    [switch]$JsonSummary,
    [string]$SummaryPath = ""
)

$ErrorActionPreference = "Stop"

Set-Location $PSScriptRoot
Set-Location ..

$projectRoot = (Get-Location).Path
$resolvedOutput = Join-Path $projectRoot $OutputDir
$env:PYTHONPATH = Join-Path $projectRoot "src"

New-Item -ItemType Directory -Path $resolvedOutput -Force | Out-Null

$startTime = Get-Date
$script:Report = [ordered]@{
    started = $startTime.ToString("o")
    project_root = $projectRoot
    output_dir = $resolvedOutput
    checks_skipped = [bool]$SkipChecks
    status = "PASS"
    steps = @()
}

function Add-Step {
    param(
        [string]$Name,
        [string]$Status,
        [string]$Message = ""
    )

    $script:Report["steps"] += [ordered]@{
        step = $Name
        status = $Status
        message = $Message
        at = (Get-Date).ToString("o")
    }
}

function Write-Step {
    param(
        [string]$Text,
        [string]$Status = "INFO"
    )

    switch ($Status) {
        "PASS" { Write-Host "  [$Status] $Text" -ForegroundColor Green }
        "WARN" { Write-Host "  [$Status] $Text" -ForegroundColor Yellow }
        "FAIL" { Write-Host "  [$Status] $Text" -ForegroundColor Red }
        default { Write-Host "  $Text" -ForegroundColor Gray }
    }
}

function Invoke-Checked {
    param(
        [string]$Name,
        [scriptblock]$Body
    )

    try {
        & $Body
        if ($LASTEXITCODE -ne 0) {
            throw "ExitCode=$LASTEXITCODE"
        }
        Add-Step -Name $Name -Status "PASS" -Message "ExitCode=0"
        Write-Step -Text "$Name (ok)" -Status "PASS"
    } catch {
        Add-Step -Name $Name -Status "FAIL" -Message $_.Exception.Message
        Write-Step -Text "$Name (failed)" -Status "FAIL"
        throw
    }
}

function Get-PythonExecutable {
    $candidates = @(
        (Join-Path $projectRoot ".venv\Scripts\python.exe"),
        (Join-Path $projectRoot ".venv\bin\python")
    )
    foreach ($candidate in $candidates) {
        if (Test-Path $candidate) {
            return $candidate
        }
    }

    $pythonCommand = Get-Command python -ErrorAction SilentlyContinue
    if ($null -ne $pythonCommand) {
        return $pythonCommand.Source
    }

    throw "Python wurde nicht gefunden. Bitte Python/venv im Pfad aktivieren."
}

$python = Get-PythonExecutable

function Invoke-Ontseq {
    param(
        [Parameter(Mandatory, Position = 0, ValueFromRemainingArguments = $true)]
        [string[]]$Arguments
    )

    if ($Arguments.Count -eq 0) {
        throw "Invoke-Ontseq benötigt mindestens ein Kommando."
    }

    $fullArgs = @($Arguments)
    if (Get-Command ontseq -ErrorAction SilentlyContinue) {
        & ontseq @fullArgs
    } else {
        & $python -m ontseq_platform.entrypoint @fullArgs
    }

    if ($LASTEXITCODE -ne 0) {
        throw "Befehl fehlgeschlagen: ontseq $($fullArgs -join ' ')"
    }
}

$failed = $false
$failMessage = ""

try {
    if (-not $NoInstall) {
        Write-Host "1/5 Installation: Editable Install inkl. dev/workflow"
        Invoke-Checked "pip install -e ." { & $python -m pip install -e ".[dev,workflow]"}
        if ($LASTEXITCODE -ne 0) { throw "Installation fehlgeschlagen." }
    } else {
        Write-Step "1/5 Installation: übersprungen (-NoInstall)" "WARN"
        Add-Step -Name "Installation" -Status "SKIP" -Message "Übersprungen via -NoInstall"
    }

    if (-not $SkipChecks) {
        Write-Host "2/5 Qualitätsprüfungen"
        Invoke-Checked "check_repository_safety.py" { & $python scripts/check_repository_safety.py }
        Invoke-Checked "check_version_consistency.py" { & $python scripts/check_version_consistency.py }
        Invoke-Checked "ruff check ." { & $python -m ruff check . }
        Invoke-Checked "ruff format --check ." { & $python -m ruff format --check . }
        Invoke-Checked "mypy src" { & $python -m mypy src }
        Invoke-Checked "unittest discover" { & $python -m unittest discover -s tests -v }
    } else {
        Write-Step "2/5 Qualitätsprüfungen: übersprungen (-SkipChecks)" "WARN"
        Add-Step -Name "Safety checks" -Status "SKIP" -Message "Übersprungen via -SkipChecks"
    }

    if (-not $SkipChecks) {
        Write-Host "3/5 Runtime Health-Check (doctor)"
        Invoke-Checked "ontseq doctor --strict --output-dir" { Invoke-Ontseq doctor --strict --output-dir $resolvedOutput }
    } else {
        Write-Step "3/5 Runtime Health-Check: übersprungen (-SkipChecks)" "WARN"
        Add-Step -Name "Runtime doctor" -Status "SKIP" -Message "Übersprungen via -SkipChecks"
    }

    Write-Host "4/5 Demo-Ausführung (synthetisch, ohne klinische Daten)"
    Invoke-Checked "ontseq demo" {
        Invoke-Ontseq demo --output-dir (Join-Path $resolvedOutput "demo")
    }

    Write-Host "5/5 System-Smoke (nur falls Toolchain vollständig installiert)"
    Write-Host "Übersprungen: optionaler Schritt. Starte bei Bedarf manuell:"
    Write-Host "  python -m ontseq_platform.entrypoint system-smoke --output-dir $([System.IO.Path]::Combine($resolvedOutput, 'system-smoke'))"
    Add-Step -Name "system-smoke" -Status "SKIP" -Message "optional"
} catch {
    $failed = $true
    $failMessage = $_.Exception.Message
    $script:Report["status"] = "FAIL"
    Add-Step -Name "pipeline" -Status "FAIL" -Message $failMessage
} finally {
    $script:Report["finished"] = (Get-Date).ToString("o")
    $script:Report["duration_seconds"] = [math]::Round(((Get-Date) - [DateTime]$startTime).TotalSeconds, 2)

    if ($JsonSummary -or $SummaryPath) {
        if ($SummaryPath -eq "") {
            $SummaryPath = Join-Path $resolvedOutput "run_test_suite_summary.json"
        }
        $json = ($script:Report | ConvertTo-Json -Depth 8)
        $summaryDir = Split-Path $SummaryPath
        if ([string]::IsNullOrWhiteSpace($summaryDir)) {
            $summaryDir = "."
        }
        [void](New-Item -ItemType Directory -Path $summaryDir -Force)
        Set-Content -Path $SummaryPath -Value $json -Encoding UTF8
        Write-Step "JSON summary written: $SummaryPath" "PASS"
        if ($JsonSummary) {
            Write-Output $json
        }
    }
}

if ($failed) {
    Write-Host "Abbruch: $failMessage" -ForegroundColor Red
    throw $failMessage
}

Write-Host "Fertig. Prüfe die erzeugten Dateien unter:"
Write-Host "  $resolvedOutput"
