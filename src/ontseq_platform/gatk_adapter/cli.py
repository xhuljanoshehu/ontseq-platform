"""Standalone research CLI. No cloud transfers, automatic execution, or clinical release."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from .contracts import BamSample, LockedFile, Mutect2Config, ReferenceBundle, RuntimeLock, VcfResource
from .core import build_plan, config_sha256, preflight, run_mutect2
from .vcf import sha256_file


def _lock(path: Path) -> LockedFile:
    path = path.expanduser().resolve(strict=True)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Cannot lock an empty/non-file input: {path.name}')
    return LockedFile(path=path, sha256=sha256_file(path))


def _sample(bam: Path, sample_id: str, reference_sha256: str) -> BamSample:
    bam = bam.expanduser().resolve(strict=True)
    indices = [p for p in (Path(str(bam) + '.bai'), bam.with_suffix('.bai')) if p.is_file()]
    if len(indices) != 1:
        raise ValueError('Exactly one conventional BAI sidecar must accompany each BAM')
    return BamSample(sample_id=sample_id, bam=_lock(bam), bai=_lock(indices[0]),
                     reference_sha256=reference_sha256)


def _resource(path: Path, name: str, reference_digest: str,
              assay_id: str | None = None) -> VcfResource:
    path = path.expanduser().resolve(strict=True)
    return VcfResource(
        name=name, vcf=_lock(path), index=_lock(Path(str(path) + '.tbi')),
        reference_sha256=reference_digest, kind='ont_pon' if assay_id else 'population',
        assay_id=assay_id,
    )


def lock_config(
    *, run_id: str, assay_id: str, reference_id: str, genome_build: str, reference: Path,
    tumor_bam: Path, tumor_sample: str, germline: Path, gatk_jar: Path,
    normal_bam: Path | None = None, normal_sample: str | None = None,
    pon: Path | None = None, contamination_sites: Path | None = None,
    intervals: Path | None = None,
) -> Mutect2Config:
    """Fingerprint local files; reference/assay labels are operator attestations.

    This does not inspect biological validity or prove resource provenance. It deliberately
    leaves both execution opt-ins false and does not download reference resources.
    """
    if (normal_bam is None) != (normal_sample is None):
        raise ValueError('Both normal BAM and normal sample name must be provided together')
    reference = reference.expanduser().resolve(strict=True)
    fasta = _lock(reference)
    ref = ReferenceBundle(
        reference_id=reference_id, genome_build=genome_build, fasta=fasta,
        fai=_lock(Path(str(reference) + '.fai')),
        dictionary=_lock(reference.with_suffix('.dict')),
    )
    config = Mutect2Config(
        run_id=run_id, assay_id=assay_id,
        mode='tumor_normal' if normal_bam else 'tumor_only', reference=ref,
        tumor=_sample(tumor_bam, tumor_sample, fasta.sha256),
        normal=_sample(normal_bam, normal_sample, fasta.sha256) if normal_bam else None,
        germline_resource=_resource(germline, 'germline-af', fasta.sha256),
        panel_of_normals=_resource(pon, 'ont-pon', fasta.sha256, assay_id) if pon else None,
        estimate_contamination=contamination_sites is not None,
        contamination_sites=_resource(contamination_sites, 'population-sites', fasta.sha256)
        if contamination_sites else None,
        intervals=_lock(intervals) if intervals else None,
        runtime=RuntimeLock(gatk_jar=_lock(gatk_jar)),
    )
    preflight(config)
    return config


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Experimental local ONT/Mutect2 adapter. Research Use Only. No clinical release.'
    )
    actions = parser.add_subparsers(dest='action', required=True)
    locking = actions.add_parser('lock', help='Hash local inputs and create a disabled configuration')
    for name in ('run-id', 'assay-id', 'reference-id', 'tumor-sample'):
        locking.add_argument('--' + name, required=True)
    locking.add_argument('--genome-build', choices=['GRCh37', 'GRCh38', 'synthetic'], required=True)
    for name in ('reference', 'tumor-bam', 'germline', 'gatk-jar', 'output'):
        locking.add_argument('--' + name, type=Path, required=True)
    for name in ('normal-bam', 'pon', 'contamination-sites', 'intervals'):
        locking.add_argument('--' + name, type=Path)
    locking.add_argument('--normal-sample')
    for action in ('plan', 'run'):
        sub = actions.add_parser(action)
        sub.add_argument('--config', type=Path, required=True)
        sub.add_argument('--output-dir', type=Path, required=True)
        if action == 'run':
            sub.add_argument('--enable', action='store_true')
            sub.add_argument('--allow-experimental-ont', action='store_true')
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.action == 'lock':
            arguments = vars(args).copy()
            arguments.pop('action')
            output = arguments.pop('output').expanduser().absolute()
            if output.exists():
                raise FileExistsError('Refusing to overwrite an existing configuration')
            cfg = lock_config(**arguments)
            descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            with os.fdopen(descriptor, 'w', encoding='utf-8') as handle:
                handle.write(cfg.model_dump_json(indent=2) + '\n')
            print(json.dumps({'status': 'NOT_RUN', 'config_sha256': config_sha256(cfg),
                              'message': 'Local hashes recorded; execution remains disabled.'}))
            return 0
        cfg = Mutect2Config.model_validate_json(args.config.read_text(encoding='utf-8'))
        if args.action == 'plan':
            steps = build_plan(cfg, args.output_dir)
            print(json.dumps({
                'status': 'NOT_RUN', 'research_only': True, 'config_sha256': config_sha256(cfg),
                'steps': [{'name': step.name, 'argv': step.argv,
                           'outputs': [str(p) for p in step.outputs]} for step in steps],
                'notice': 'Planning only; no native tools, integrity checks or biological analysis ran.',
            }, indent=2))
            return 0
        options = cfg.model_dump(mode='json')
        if args.enable:
            options['enabled'] = True
        if args.allow_experimental_ont:
            options['allow_experimental_ont'] = True
        cfg = Mutect2Config.model_validate(options)
        result = run_mutect2(cfg, args.output_dir)
        print(json.dumps({
            'status': result.status, 'run_id': result.run_id,
            'candidate_count': len(result.candidates),
            'clinically_reportable': False, 'messages': result.messages,
        }, indent=2))
        return {'COMPLETED': 0, 'FAILED': 1, 'NOT_RUN': 2}[result.status]
    except (ValueError, OSError, KeyError) as exc:
        print(json.dumps({'status': 'FAILED', 'error': str(exc),
                          'clinically_reportable': False}), file=sys.stderr)
        return 1
