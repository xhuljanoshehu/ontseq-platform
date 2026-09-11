from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "desktop-ci.yml"


def test_desktop_bundle_builds_the_runtime_installer_contract() -> None:
    """The published bundle contains every file RuntimePackage.VerifyAsync requires."""
    workflow = WORKFLOW.read_text(encoding="utf-8")

    required_build_contract = (
        "ontseq-linux-runtime.tar.gz",
        "ontseq_platform-${{ env.ONTSEQ_VERSION }}-py3-none-any.whl",
        "SHA256SUMS",
    )
    for marker in required_build_contract:
        assert marker in workflow, (
            f"desktop-ci.yml does not package required runtime input: {marker}"
        )

    required_bundle_checks = (
        '$wheel = "desktop/publish/runtime/'
        'ontseq_platform-${{ env.ONTSEQ_VERSION }}-py3-none-any.whl"',
        '$checksums = "desktop/publish/runtime/SHA256SUMS"',
    )
    for marker in required_bundle_checks:
        assert marker in workflow, (
            f"desktop-ci.yml does not verify required bundle input: {marker}"
        )
