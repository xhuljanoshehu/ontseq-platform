from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.service.guard import (
    GuardError,
    host_is_loopback,
    new_token,
    origin_is_loopback,
    resolve_bam_index,
    resolve_within,
    token_matches,
    windows_to_wsl,
    wsl_to_windows,
)


class TokenTests(unittest.TestCase):
    def test_two_tokens_differ(self) -> None:
        self.assertNotEqual(new_token(), new_token())

    def test_a_token_matches_itself(self) -> None:
        token = new_token()
        self.assertTrue(token_matches(token, token))

    def test_a_wrong_token_is_refused(self) -> None:
        self.assertFalse(token_matches("wrong", new_token()))

    def test_a_missing_token_is_refused_rather_than_treated_as_empty(self) -> None:
        """A request carrying no token must not match a service that has one."""
        self.assertFalse(token_matches(None, new_token()))
        self.assertFalse(token_matches("", new_token()))

    def test_a_prefix_of_the_token_is_refused(self) -> None:
        token = new_token()
        self.assertFalse(token_matches(token[:-1], token))


class HostTests(unittest.TestCase):
    """The DNS-rebinding check. The socket is local in that attack; the header is not."""

    def test_loopback_names_are_accepted(self) -> None:
        for host in ("127.0.0.1:8765", "localhost:8765", "[::1]:8765"):
            self.assertTrue(host_is_loopback(host, port=8765), host)

    def test_a_foreign_name_pointed_at_loopback_is_refused(self) -> None:
        self.assertFalse(host_is_loopback("evil.example.com:8765", port=8765))

    def test_the_wrong_port_is_refused(self) -> None:
        self.assertFalse(host_is_loopback("127.0.0.1:9999", port=8765))

    def test_a_missing_host_is_refused(self) -> None:
        self.assertFalse(host_is_loopback(None, port=8765))

    def test_a_bare_loopback_name_without_a_port_is_accepted(self) -> None:
        self.assertTrue(host_is_loopback("localhost", port=8765))


class OriginTests(unittest.TestCase):
    def test_an_absent_origin_is_allowed_because_browsers_omit_it(self) -> None:
        self.assertTrue(origin_is_loopback(None, port=8765))

    def test_our_own_origin_is_allowed(self) -> None:
        self.assertTrue(origin_is_loopback("http://127.0.0.1:8765", port=8765))

    def test_another_tab_is_refused(self) -> None:
        self.assertFalse(origin_is_loopback("https://example.com", port=8765))

    def test_a_null_origin_is_refused(self) -> None:
        """A file:// page sends null. It is not this page, so it does not get in."""
        self.assertFalse(origin_is_loopback("null", port=8765))

    def test_nonsense_is_refused_rather_than_parsed_generously(self) -> None:
        for origin in ("", "127.0.0.1:8765", "javascript:alert(1)"):
            self.assertFalse(origin_is_loopback(origin, port=8765), origin)


class RootBoundaryTests(unittest.TestCase):
    def test_a_path_inside_an_allowed_root_resolves(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "runs" / "sample.bam"
            target.parent.mkdir()
            target.write_bytes(b"x")
            self.assertEqual(resolve_within(target, [root]), target.resolve())

    def test_a_path_outside_every_root_is_refused(self) -> None:
        with (
            tempfile.TemporaryDirectory() as a,
            tempfile.TemporaryDirectory() as b,
            self.assertRaises(GuardError),
        ):
            resolve_within(Path(b) / "secret", [Path(a)])

    def test_dot_dot_cannot_climb_out(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "allowed"
            root.mkdir()
            with self.assertRaises(GuardError):
                resolve_within(root / ".." / ".." / "etc" / "passwd", [root])

    def test_a_symlink_pointing_out_of_the_root_is_refused(self) -> None:
        """The obvious way past a prefix check, so resolution happens before comparison."""
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root, outside = base / "allowed", base / "outside"
            root.mkdir()
            outside.mkdir()
            (outside / "secret.bam").write_bytes(b"x")
            (root / "link.bam").symlink_to(outside / "secret.bam")
            with self.assertRaises(GuardError):
                resolve_within(root / "link.bam", [root])

    def test_no_roots_means_nothing_is_readable(self) -> None:
        """A service that was never told what it may read must not read everything."""
        with self.assertRaises(GuardError):
            resolve_within("/etc/passwd", [])

    def test_the_root_itself_is_inside_itself(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.assertEqual(resolve_within(root, [root]), root.resolve())

    def test_a_sibling_with_the_same_prefix_is_not_inside(self) -> None:
        """`/data/runs-old` must not pass because `/data/runs` is allowed."""
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            allowed, sibling = base / "runs", base / "runs-old"
            allowed.mkdir()
            sibling.mkdir()
            with self.assertRaises(GuardError):
                resolve_within(sibling, [allowed])


class BamIndexResolutionTests(unittest.TestCase):
    def test_bam_dot_bai_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            index = Path(f"{bam}.bai")
            bam.write_bytes(b"BAM")
            index.write_bytes(b"BAI")
            self.assertEqual(resolve_bam_index(bam), index)

    def test_short_bai_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            index = bam.with_suffix(".bai")
            bam.write_bytes(b"BAM")
            index.write_bytes(b"BAI")
            self.assertEqual(resolve_bam_index(bam), index)

    def test_bam_dot_bai_has_deterministic_precedence_when_both_exist(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            preferred = Path(f"{bam}.bai")
            alternative = bam.with_suffix(".bai")
            for path in (bam, preferred, alternative):
                path.write_bytes(b"x")
            self.assertEqual(resolve_bam_index(bam), preferred)

    def test_missing_index_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            bam.write_bytes(b"BAM")
            with self.assertRaises(GuardError) as caught:
                resolve_bam_index(bam)
            self.assertIn("sample.bam.bai", str(caught.exception))
            self.assertIn("sample.bai", str(caught.exception))

    def test_unrelated_bai_is_not_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "sample.bam"
            bam.write_bytes(b"BAM")
            (root / "other.bai").write_bytes(b"BAI")
            with self.assertRaises(GuardError):
                resolve_bam_index(bam)


class PathTranslationTests(unittest.TestCase):
    def test_a_drive_letter_becomes_a_mount_point(self) -> None:
        self.assertEqual(windows_to_wsl(r"P:\Lab_FG06\run.bam"), "/mnt/p/Lab_FG06/run.bam")

    def test_forward_slashes_are_accepted_too(self) -> None:
        self.assertEqual(windows_to_wsl("C:/Users/x/a.bam"), "/mnt/c/Users/x/a.bam")

    def test_a_posix_path_passes_through_unchanged(self) -> None:
        self.assertEqual(windows_to_wsl("/mnt/p/run.bam"), "/mnt/p/run.bam")

    def test_surrounding_quotes_from_a_paste_are_stripped(self) -> None:
        self.assertEqual(windows_to_wsl('"P:\\a\\b.bam"'), "/mnt/p/a/b.bam")

    def test_a_unc_path_is_refused_with_the_reason(self) -> None:
        """WSL does not map these on its own; inventing a mount point invents a path."""
        with self.assertRaises(GuardError) as caught:
            windows_to_wsl(r"\\fileserver\Lab_FG06\run.bam")
        self.assertIn("mount", str(caught.exception))

    def test_the_round_trip_returns_the_original(self) -> None:
        self.assertEqual(wsl_to_windows(windows_to_wsl(r"P:\Lab\run.bam")), r"P:\Lab\run.bam")

    def test_a_path_that_is_not_a_windows_mount_is_left_alone(self) -> None:
        for path in ("/home/grid/run.bam", "/mnt/data/run.bam"):
            self.assertEqual(wsl_to_windows(path), path)


if __name__ == "__main__":
    unittest.main()
