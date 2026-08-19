from __future__ import annotations


def main() -> None:
    try:
        from .desktop.app import run
    except ImportError as exc:
        raise SystemExit(
            'ONTSeq Desktop requires the desktop extra. Install with: pip install -e ".[desktop]"'
        ) from exc
    raise SystemExit(run())


if __name__ == "__main__":
    main()
