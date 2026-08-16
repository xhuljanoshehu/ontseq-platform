from __future__ import annotations

import argparse
import json
from pathlib import Path

from ontseq_platform.cnv.cytobands import CytobandTable
from ontseq_platform.cnv.models import (
    CnvBenchmarkCase,
    CnvCallSet,
    CnvEvaluationReport,
    CnvTruthSet,
)
from ontseq_platform.cnv.strata import CnvAggregateReport
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    BenchmarkCase,
    BenchmarkReport,
    CraminoQCReport,
    LocalSmokeReport,
    PipelineResult,
    ReferenceLock,
    SampleManifest,
    SnifflesCallReport,
    SnifflesPolicy,
)


def _render() -> dict[Path, str]:
    return {
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
        Path("schemas/cnv-truth-set.schema.json"): json.dumps(
            CnvTruthSet.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/cnv-call-set.schema.json"): json.dumps(
            CnvCallSet.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/cnv-benchmark-case.schema.json"): json.dumps(
            CnvBenchmarkCase.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/cnv-evaluation-report.schema.json"): json.dumps(
            CnvEvaluationReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/cnv-aggregate-report.schema.json"): json.dumps(
            CnvAggregateReport.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/cytoband-table.schema.json"): json.dumps(
            CytobandTable.model_json_schema(), indent=2, sort_keys=True
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
            path.write_text(content, encoding="utf-8")
            print(path)
    if mismatches:
        raise SystemExit(f"Schemas are stale: {', '.join(mismatches)}")


if __name__ == "__main__":
    main()
