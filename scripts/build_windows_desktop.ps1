param(
    [switch]$SkipTests
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $RepoRoot

Write-Host "Building ONTSeq Desktop from $RepoRoot"
python -m pip install --upgrade pip
python -m pip install -e ".[desktop,dev]"

if (-not $SkipTests) {
    python -m unittest discover -s tests -v
    ruff check .
    ruff format --check .
}

Remove-Item -Recurse -Force build -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force dist -ErrorAction SilentlyContinue
Remove-Item -Force ONTSeqDesktop.spec -ErrorAction SilentlyContinue

python -m PyInstaller `
    --noconfirm `
    --clean `
    --onefile `
    --windowed `
    --name ONTSeqDesktop `
    --paths src `
    packaging/desktop_launcher.py

$Exe = Join-Path $RepoRoot "dist\ONTSeqDesktop.exe"
if (-not (Test-Path $Exe)) {
    throw "Expected executable was not created: $Exe"
}

Write-Host "Created: $Exe"
