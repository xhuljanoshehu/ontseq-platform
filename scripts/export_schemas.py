from __future__ import annotations

import argparse
import json
from pathlib import Path

from ontseq_platform.dilution import (
    DilutionPolicy,
    DilutionSeriesPlan,
    DilutionSeriesReport,
    LodPolicy,
    LodReport,
)
from ontseq_platform.methylation import MethylationPolicy, MethylationReport
from ontseq_platform.methylation_holdout import (
    PairedHoldoutCohort,
    PairedHoldoutEvidence,
    PairedHoldoutRecoveryReport,
    PairedHoldoutRegistration,
    PairedHoldoutSplitPolicy,
)
from ontseq_platform.methylation_holdout_cli import PairedHoldoutReadiness
from ontseq_platform.methylation_holdout_runner import PairedHoldoutLocalInputs
from ontseq_platform.methylation_mixture import (
    MethylationMixturePolicy,
    MethylationMixtureReport,
    NanopolishSourceMetadata,
)
from ontseq_platform.methylation_validation import (
    MethylationValidationEvidence,
    MethylationValidationMatrix,
    MethylationValidationRegistration,
    MethylationValidationReport,
    ValidationCohort,
)
from ontseq_platform.methylation_validation_runner import ValidationLocalInputs, ValidationReadiness
from ontseq_platform.modbam import (
    ModbamAdapterPolicy,
    ModbamSourceMetadata,
    ModbamSourceSummary,
    ModkitExtractCommand,
)
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AmlKnowledgeLock,
    BenchmarkCase,
    BenchmarkReport,
    CraminoQCReport,
    CuteSvCallReport,
    CuteSvPolicy,
    IntervalResourceLock,
    LocalSmokeReport,
    PipelineResult,
    ReferenceLock,
    SampleManifest,
    SnifflesCallReport,
    SnifflesPolicy,
    SvConsensusPolicy,
    SvConsensusReport,
    SvEvidencePolicy,
)


def _render() -> dict[Path, str]:
    return {
        **{
            Path(f"schemas/{name}.schema.json"): json.dumps(
                model.model_json_schema(), indent=2, sort_keys=True
            )
            + "\n"
            for name, model in {
                "modbam-adapter-policy": ModbamAdapterPolicy,
                "modbam-source-metadata": ModbamSourceMetadata,
                "modbam-source-summary": ModbamSourceSummary,
                "modkit-extract-command": ModkitExtractCommand,
                "methylation-validation-matrix": MethylationValidationMatrix,
                "methylation-validation-cohort": ValidationCohort,
                "methylation-validation-registration": MethylationValidationRegistration,
                "methylation-validation-evidence": MethylationValidationEvidence,
                "methylation-validation-report": MethylationValidationReport,
                "methylation-validation-local-inputs": ValidationLocalInputs,
                "methylation-validation-readiness": ValidationReadiness,
                "methylation-holdout-cohort": PairedHoldoutCohort,
                "methylation-holdout-evidence": PairedHoldoutEvidence,
                "methylation-holdout-report": PairedHoldoutRecoveryReport,
                "methylation-holdout-registration": PairedHoldoutRegistration,
                "methylation-holdout-split-policy": PairedHoldoutSplitPolicy,
                "methylation-holdout-readiness": PairedHoldoutReadiness,
                "methylation-holdout-local-inputs": PairedHoldoutLocalInputs,
            }.items()
        },
        Path("schemas/sample-manifest.schema.json"): json.dumps(
            SampleManifest.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/pipeline-result.schema.json"): json.dumps(
            PipelineResult.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/reference-lock.schema.json"): json.dumps(
            ReferenceLock.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/aligned-bam-intake.schema.json"): json.dumps(
            AlignedBamIntakeReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/cramino-qc.schema.json"): json.dumps(
            CraminoQCReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/sniffles-policy.schema.json"): json.dumps(
            SnifflesPolicy.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/sniffles-call.schema.json"): json.dumps(
            SnifflesCallReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/cutesv-policy.schema.json"): json.dumps(
            CuteSvPolicy.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/cutesv-call.schema.json"): json.dumps(
            CuteSvCallReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/sv-consensus-policy.schema.json"): json.dumps(
            SvConsensusPolicy.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/sv-consensus-report.schema.json"): json.dumps(
            SvConsensusReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/sv-evidence-policy.schema.json"): json.dumps(
            SvEvidencePolicy.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/interval-resource-lock.schema.json"): json.dumps(
            IntervalResourceLock.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/aml-knowledge-lock.schema.json"): json.dumps(
            AmlKnowledgeLock.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/local-smoke.schema.json"): json.dumps(
            LocalSmokeReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/benchmark-case.schema.json"): json.dumps(
            BenchmarkCase.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/benchmark-report.schema.json"): json.dumps(
            BenchmarkReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/methylation-policy.schema.json"): json.dumps(
            MethylationPolicy.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/methylation-report.schema.json"): json.dumps(
            MethylationReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/methylation-mixture-policy.schema.json"): json.dumps(
            MethylationMixturePolicy.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/methylation-mixture-report.schema.json"): json.dumps(
            MethylationMixtureReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/nanopolish-source-metadata.schema.json"): json.dumps(
            NanopolishSourceMetadata.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/dilution-policy.schema.json"): json.dumps(
            DilutionPolicy.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/dilution-series-plan.schema.json"): json.dumps(
            DilutionSeriesPlan.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/dilution-series-report.schema.json"): json.dumps(
            DilutionSeriesReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/lod-policy.schema.json"): json.dumps(
            LodPolicy.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/lod-report.schema.json"): json.dumps(
            LodReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    mismatches: list[str] = []
    for path, content in _render().items():
        if args.check:
            if not path.exists() or path.read_text(encoding="utf-8") != content:
                mismatches.append(str(path))
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8", newline="\n")
            print(path)
    if mismatches:
        raise SystemExit(f"Schemas are stale: {', '.join(mismatches)}")


if __name__ == "__main__":
    main()
