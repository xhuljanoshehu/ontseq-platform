@echo off
rem ============================================================================
rem  ONTSeq — Start von Windows aus.
rem
rem  Auf den Desktop legen und doppelklicken. Der Rest passiert von selbst:
rem  der Dienst startet in WSL, der Browser geht auf.
rem
rem  EINMALIG ANPASSEN — die drei Zeilen unten, dann nie wieder.
rem  Pfade in WSL-Schreibweise: aus P:\NANOPORE wird /mnt/p/NANOPORE.
rem ============================================================================

set "DISTRO=Ubuntu"
set "REFLOCK=/mnt/p/NANOPORE/hg38.reference-lock.json"
set "DATEN=/mnt/p/NANOPORE"
set "AUSGABE=/mnt/p/AUSWERTUNG"
set "PORT=8765"

rem ---------------------------------------------------------------------------
title ONTSeq
echo.
echo   ONTSeq wird gestartet ...
echo.

rem Laeuft der Dienst schon? Dann nur das Fenster oeffnen. Ein zweiter Start
rem wuerde am belegten Port scheitern und wie ein Fehler aussehen, obwohl alles
rem in Ordnung ist.
powershell -NoProfile -Command ^
  "try { $null = Invoke-WebRequest -Uri 'http://127.0.0.1:%PORT%/' -TimeoutSec 2 -UseBasicParsing; exit 0 } catch { exit 1 }" >nul 2>&1
if %ERRORLEVEL% EQU 0 (
  echo   Der Dienst laeuft bereits.
  start "" "http://127.0.0.1:%PORT%/"
  timeout /t 2 >nul
  exit /b 0
)

rem Ist WSL ueberhaupt da? Sonst kommt sonst eine Fehlermeldung, die niemandem
rem sagt, was zu tun ist.
wsl -d %DISTRO% -- true >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
  echo   FEHLER: WSL-Distribution "%DISTRO%" ist nicht erreichbar.
  echo.
  echo   Pruefen Sie, ob WSL installiert ist ^(wsl --list --verbose^) und ob der
  echo   Name oben in dieser Datei stimmt.
  echo.
  pause
  exit /b 1
)

start "" "http://127.0.0.1:%PORT%/"

echo   Dieses Fenster bitte offen lassen — es haelt den Dienst am Laufen.
echo   Zum Beenden dieses Fenster schliessen.
echo.

wsl -d %DISTRO% -- bash -lc "ontseq serve --reference-lock '%REFLOCK%' --allow-root '%DATEN%' --output-dir '%AUSGABE%' --port %PORT% --no-browser"

echo.
echo   Der Dienst wurde beendet.
pause
