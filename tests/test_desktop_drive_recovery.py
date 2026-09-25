"""Exercise the embedded recovery helper without changing any real mounts."""

import errno
import importlib.util
import subprocess
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[1] / "desktop/ONTSeq.Desktop/wsl_drive_recovery.py"


@pytest.fixture
def recovery():
    assert SOURCE.is_file(), "Desktop recovery helper is not implemented"
    spec = importlib.util.spec_from_file_location("drive_recovery", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mount_line(source="E:\\134", kind="9p", options="rw,aname=drvfs;path=E:\\134"):
    return f"22 1 0:2 / /mnt/e rw - {kind} {source} {options}\n"


def configure(monkeypatch, recovery, *, err=errno.ENODEV, mount=None, busy=False):
    calls = []
    state = {"stale": err}

    def stat(path):
        if state["stale"]:
            raise OSError(state["stale"], "synthetic drive error")
        return None

    def run(argv, **kwargs):
        calls.append(argv)
        if argv[0] == "/bin/umount" and busy:
            raise subprocess.CalledProcessError(32, argv, stderr="target is busy")
        if argv[0] == "/bin/mount":
            state["stale"] = None
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(recovery, "stat_mount", stat)
    monkeypatch.setattr(
        recovery, "read_mounts", lambda: mount if mount is not None else mount_line()
    )
    monkeypatch.setattr(recovery.subprocess, "run", run)
    return calls


def test_stale_drive_is_detected_without_changes_then_repaired(monkeypatch, recovery):
    calls = configure(monkeypatch, recovery)
    assert recovery.recover("E", 1000, 1000, repair=False)["state"] == "stale"
    assert calls == []
    assert recovery.recover("E", 1000, 1000, repair=True)["state"] == "repaired"
    assert calls == [
        ["/bin/umount", "/mnt/e"],
        ["/bin/mount", "-t", "drvfs", "E:", "/mnt/e", "-o", "uid=1000,gid=1000"],
    ]


@pytest.mark.parametrize("err", [None, errno.EACCES, errno.ENOENT, errno.EIO])
def test_healthy_and_unrelated_errors_never_unmount(monkeypatch, recovery, err):
    calls = configure(monkeypatch, recovery, err=err)
    if err is None:
        assert recovery.recover("E", 1000, 1000, repair=True)["state"] == "healthy"
    else:
        with pytest.raises(OSError):
            recovery.recover("E", 1000, 1000, repair=True)
    assert calls == []


@pytest.mark.parametrize("drive", ["C", "/", "E;shutdown", "EE", "../e"])
def test_system_drive_and_invalid_drive_are_refused(monkeypatch, recovery, drive):
    calls = configure(monkeypatch, recovery)
    with pytest.raises(ValueError):
        recovery.recover(drive, 1000, 1000, repair=True)
    assert calls == []


@pytest.mark.parametrize(
    "mount",
    [
        mount_line(source="F:\\134"),
        mount_line(kind="ext4"),
        "",
        mount_line(options="rw,aname=other"),
    ],
)
def test_unknown_or_mismatched_mount_is_not_replaced(monkeypatch, recovery, mount):
    calls = configure(monkeypatch, recovery, mount=mount)
    with pytest.raises(ValueError):
        recovery.recover("E", 1000, 1000, repair=True)
    assert calls == []


def test_busy_mount_is_not_forced_or_replaced(monkeypatch, recovery):
    calls = configure(monkeypatch, recovery, busy=True)
    with pytest.raises(subprocess.CalledProcessError):
        recovery.recover("E", 1000, 1000, repair=True)
    assert calls == [["/bin/umount", "/mnt/e"]]


def test_drvfs_source_without_trailing_backslash_is_supported(monkeypatch, recovery):
    calls = configure(monkeypatch, recovery, mount=mount_line(source="E:"))
    assert recovery.recover("E", 1000, 1000, repair=True)["state"] == "repaired"
    assert calls[0] == ["/bin/umount", "/mnt/e"]


def test_unmounted_directory_is_not_reported_healthy(monkeypatch, recovery):
    calls = configure(monkeypatch, recovery, err=None, mount="")
    with pytest.raises(ValueError):
        recovery.recover("E", 1000, 1000, repair=False)
    assert calls == []


@pytest.mark.parametrize(
    "options",
    [
        "rw,aname=drvfs;path=E:\\134;metadata",
        "rw,aname=drvfs;path=E:\\134;umask=077",
        "rw,aname=drvfs;path=E:\\134;fmask=011",
        "rw,aname=drvfs;path=E:\\134;dmask=77",
        "rw,aname=drvfs;path=E:\\134;case=dir",
        "rw,aname=drvfs;path=E:\\134;case=force",
    ],
)
def test_nondefault_drvfs_semantics_fail_closed_before_unmount(monkeypatch, recovery, options):
    calls = configure(monkeypatch, recovery, mount=mount_line(options=options))
    assert recovery.recover("E", 1000, 1000, repair=False)["state"] == "stale"
    with pytest.raises(ValueError, match="nicht sicher rekonstruiert"):
        recovery.recover("E", 1000, 1000, repair=True)
    assert calls == []


@pytest.mark.parametrize(
    "options",
    [
        "rw,aname=drvfs;path=E:\\134;uid=0;gid=1000",
        "rw,aname=drvfs;path=E:\\134;uid=1000;gid=0",
    ],
)
def test_mismatched_drvfs_identity_fails_closed_before_unmount(monkeypatch, recovery, options):
    calls = configure(monkeypatch, recovery, mount=mount_line(options=options))
    with pytest.raises(ValueError, match="WSL-Benutzer"):
        recovery.recover("E", 1000, 1000, repair=True)
    assert calls == []


def test_equivalent_explicit_drvfs_identity_and_default_case_remain_repairable(
    monkeypatch, recovery
):
    calls = configure(
        monkeypatch,
        recovery,
        mount=mount_line(
            options="rw,aname=drvfs;path=E:\\134;uid=1000;gid=1000;umask=022;fmask=000;dmask=000;case=off"
        ),
    )
    assert recovery.recover("E", 1000, 1000, repair=True)["state"] == "repaired"
    assert calls == [
        ["/bin/umount", "/mnt/e"],
        ["/bin/mount", "-t", "drvfs", "E:", "/mnt/e", "-o", "uid=1000,gid=1000"],
    ]
