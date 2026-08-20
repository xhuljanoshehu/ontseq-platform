"""Dorado basecalling adapter.

.. warning::

   **This adapter has never been executed against a real Dorado binary**, here or in
   continuous integration. Dorado needs a GPU, downloaded model files and real POD5
   signal, none of which exist in the environments available to this repository's tests.

   The code is structured to the same adapter boundary as the verified adapters, and it
   fails closed on every gate it can check locally, but its behaviour against a real
   Dorado installation is an *assumption*. Every artifact it produces is marked
   ``unverified_adapter`` so that a downstream reader is never left guessing which parts
   of a run rest on tested code. The first real execution should be treated as a test of
   this adapter, not as a production run.

Why the gates below exist
-------------------------

Basecalling is the one stage where a silent misconfiguration is unrecoverable: the reads
are wrong, everything downstream is confidently wrong, and nothing later in the pipeline
can detect it. So the adapter refuses to proceed on a version mismatch, on an
unidentifiable model, and on a model directory whose checksum does not match the lock.

Modified-base calling is requested explicitly rather than left to a default, because a run
basecalled without ``MM``/``ML`` tags cannot be used for methylation later and the loss is
invisible until someone tries.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .execution import StreamingCommandRunner, SubprocessRunner
from .model_lock import fingerprint as fingerprint_model
from .models import FileFingerprint, StrictModel, ToolRecord
from .reference import sha256_file

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+\.\d+)(?!\d)")

#: Marker carried by every artifact this module produces.
ADAPTER_VERIFICATION = "unverified_adapter"

UNVERIFIED_NOTICE = (
    "Produced by an adapter that has never been executed against a real Dorado binary. "
    "Treat the first real run as a test of this adapter."
)


class BasecallPolicy(StrictModel):
    """Version- and model-locked basecalling configuration."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    expected_version: str = Field(pattern=r"^\d+\.\d+\.\d+$")
    #: Dorado model name or an absolute path to a downloaded model directory.
    model: str = Field(min_length=1)
    #: SHA-256 of the model archive or directory listing, when the site records one.
    model_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    #: Modified-base models to call alongside the basecall, for example ``5mCG_5hmCG``.
    modified_bases: list[str] = Field(default_factory=list)
    device: str = Field(default="cuda:all", min_length=1)
    #: Reads below this mean quality are dropped by Dorado itself.
    minimum_qscore: int | None = Field(default=None, ge=0, le=60)
    emit_moves: bool = False
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def modified_bases_are_unique(self) -> BasecallPolicy:
        if len(self.modified_bases) != len(set(self.modified_bases)):
            raise ValueError("modified base models must be unique")
        if self.status == "validated":
            raise ValueError(
                "this adapter has not been executed against a real Dorado binary and "
                "must not be marked validated"
            )
        return self


class BasecallReport(StrictModel):
    """Normalized record of one basecalling run."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    policy: BasecallPolicy
    tool: ToolRecord
    #: Envelope-relative name of the produced unaligned BAM.
    unaligned_bam_relative_path: str = Field(min_length=1)
    unaligned_bam_fingerprint: FileFingerprint
    pod5_file_count: int = Field(ge=1)
    modified_bases_requested: list[str] = Field(default_factory=list)
    adapter_verification: Literal["unverified_adapter"] = "unverified_adapter"
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True


@dataclass(frozen=True)
class BasecallInputs:
    """Everything the basecalling stage needs from outside the envelope."""

    pod5_directory: Path


def dorado_version(text: str) -> str:
    """Parse a Dorado version from its probe output.

    Public because preflight has to reach the same answer this module will. A preflight
    that parsed versions differently from the run it precedes could clear a run that then
    fails on the version lock, which is worse than not checking at all.
    """
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    raise ValueError("could not parse a version from dorado output")


def model_signature(model: str) -> str | None:
    """Fingerprint a local model directory, or return ``None`` for a named model.

    A directory is hashed from its sorted relative file names and their individual
    checksums, so the value is stable across machines and detects a partially downloaded
    or tampered model. A bare model name cannot be fingerprinted without resolving
    Dorado's own cache, and inventing a value would be worse than admitting none.

    Delegates to :mod:`ontseq_platform.model_lock`, which is also what ``ontseq model-lock``
    calls. Two implementations of one digest is the shape where a preflight check starts
    disagreeing with the command that produced the value it compares against.
    """
    path = Path(model)
    if not path.is_dir():
        return None
    return fingerprint_model(path).signature


def build_basecaller_argv(
    *,
    dorado: str,
    policy: BasecallPolicy,
    pod5_directory: Path,
    threads: int | None = None,
) -> list[str]:
    """Build the basecalling command as an explicit argument vector.

    Dorado writes BAM to stdout, so no output path appears here; the caller streams it to
    a file through :meth:`SubprocessRunner.run_to_file`.
    """
    argv = [dorado, "basecaller", policy.model, str(pod5_directory)]
    argv.extend(["--device", policy.device])
    for modification in policy.modified_bases:
        argv.extend(["--modified-bases", modification])
    if policy.minimum_qscore is not None:
        argv.extend(["--min-qscore", str(policy.minimum_qscore)])
    if policy.emit_moves:
        argv.append("--emit-moves")
    if threads is not None:
        argv.extend(["--threads", str(threads)])
    return argv


def run_basecalling(
    inputs: BasecallInputs,
    policy: BasecallPolicy,
    *,
    sample_id: str,
    output_bam: Path,
    runner: StreamingCommandRunner | None = None,
    dorado: str = "dorado",
    threads: int | None = None,
    timeout_seconds: int = 86400,
) -> BasecallReport:
    """Basecall a POD5 directory into an unaligned BAM.

    Fails closed before doing any work when the POD5 directory is empty, the Dorado
    version does not match the lock, or a locked model checksum does not match what is on
    disk.
    """
    if not inputs.pod5_directory.is_dir():
        raise ValueError("POD5 input directory is missing or unreadable")
    pod5_files = sorted(inputs.pod5_directory.rglob("*.pod5"))
    if not pod5_files:
        raise ValueError("POD5 input directory contains no .pod5 files")
    if output_bam.exists():
        raise ValueError("refusing to overwrite an existing basecalled BAM")

    command_runner = runner or SubprocessRunner()
    probe = command_runner.run([dorado, "--version"], timeout_seconds=120)
    if probe.returncode != 0:
        raise ValueError(f"dorado version probe returned exit code {probe.returncode}")
    version = dorado_version(f"{probe.stdout}\n{probe.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"dorado version {version!r} does not match the policy lock {policy.expected_version!r}"
        )

    observed_model = model_signature(policy.model)
    if policy.model_sha256 is not None:
        if observed_model is None:
            raise ValueError(
                "the policy locks a model checksum but the model is not a local "
                "directory that can be fingerprinted"
            )
        if observed_model != policy.model_sha256:
            raise ValueError("the basecalling model on disk does not match the policy lock")

    result = command_runner.run_to_file(
        build_basecaller_argv(
            dorado=dorado,
            policy=policy,
            pod5_directory=inputs.pod5_directory,
            threads=threads,
        ),
        output_bam,
        timeout_seconds=timeout_seconds,
    )
    if result.returncode != 0:
        detail = (result.stderr or "").strip().splitlines()
        tail = detail[-1] if detail else "no diagnostic output"
        raise ValueError(f"dorado basecaller failed with exit code {result.returncode}: {tail}")
    if not output_bam.is_file() or output_bam.stat().st_size == 0:
        raise ValueError("dorado basecaller reported success but produced no output")

    warnings = [UNVERIFIED_NOTICE]
    if not policy.modified_bases:
        warnings.append(
            "No modified-base model was requested. The resulting reads carry no MM/ML "
            "tags and cannot support a later methylation analysis."
        )
    if observed_model is None and policy.model_sha256 is None:
        warnings.append(
            "The model was given by name, so no model checksum is recorded. Provenance "
            "for this run cannot prove which model weights were used."
        )

    return BasecallReport(
        sample_id=sample_id,
        policy=policy,
        tool=ToolRecord(
            name="dorado",
            version=version,
            parameters={
                "model": policy.model,
                "model_sha256": observed_model,
                "device": policy.device,
                "modified_bases": list(policy.modified_bases),
                "minimum_qscore": policy.minimum_qscore,
                "emit_moves": policy.emit_moves,
                "adapter_verification": ADAPTER_VERIFICATION,
            },
        ),
        unaligned_bam_relative_path=output_bam.name,
        unaligned_bam_fingerprint=FileFingerprint(
            size_bytes=output_bam.stat().st_size, sha256=sha256_file(output_bam)
        ),
        pod5_file_count=len(pod5_files),
        modified_bases_requested=list(policy.modified_bases),
        warnings=warnings,
        limitations=[
            UNVERIFIED_NOTICE,
            "Basecalling accuracy depends on the model, chemistry and flow-cell version; "
            "none of that is validated here.",
            "No read-level quality gate is applied beyond Dorado's own min-qscore.",
            "Duplex basecalling, demultiplexing and adapter trimming are out of scope for "
            "this adapter.",
        ],
    )
