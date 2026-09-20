"""Synthetic contract tests; these do NOT run GATK or validate ONT performance."""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


def test_adapter_module_exists():
    assert importlib.util.find_spec('ontseq_platform.gatk_adapter.core') is not None


def locked(path: Path, text: str = 'synthetic-only\n') -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


@pytest.fixture
def config_dict(tmp_path):
    ref = locked(tmp_path / 'reference.fa', '>chrSynthetic\nACGTACGT\n')
    tumor = locked(tmp_path / 'tumor.bam')
    return {
        'run_id': 'synthetic-run-001', 'assay_id': 'synthetic-ont-assay-v1',
        'mode': 'tumor_only',
        'reference': {
            'reference_id': 'synthetic-ref-v1', 'genome_build': 'synthetic',
            'fasta': ref,
            'fai': locked(tmp_path / 'reference.fa.fai', 'chrSynthetic\t8\t14\t8\t9\n'),
            'dictionary': locked(tmp_path / 'reference.dict',
                                 '@HD\tVN:1.6\n@SQ\tSN:chrSynthetic\tLN:8\n'),
        },
        'tumor': {
            'sample_id': 'SYNTHETIC_T', 'bam': tumor,
            'bai': locked(tmp_path / 'tumor.bam.bai'),
            'reference_sha256': ref['sha256'],
        },
        'germline_resource': {
            'name': 'synthetic-germline',
            'vcf': locked(tmp_path / 'germline.vcf.gz'),
            'index': locked(tmp_path / 'germline.vcf.gz.tbi'),
            'reference_sha256': ref['sha256'], 'kind': 'population',
        },
        'runtime': {
            'gatk_jar': locked(tmp_path / 'gatk-package-4.6.2.0-local.jar'),
            'expected_version': '4.6.2.0',
        },
    }


def config(data):
    from ontseq_platform.gatk_adapter.contracts import Mutect2Config
    return Mutect2Config.model_validate(data)


def with_normal(data, tmp_path):
    data['mode'] = 'tumor_normal'
    data['normal'] = {
        'sample_id': 'SYNTHETIC_N',
        'bam': locked(tmp_path / 'normal.bam', 'different-synthetic-normal\n'),
        'bai': locked(tmp_path / 'normal.bam.bai'),
        'reference_sha256': data['reference']['fasta']['sha256'],
    }
    return data


def vcf_text(samples=('SYNTHETIC_T',), body=None):
    if body is None:
        body = 'chrSynthetic\t2\t.\tC\tT\t.\tPASS\tTLOD=8\tGT:AD:AF:DP\t0/1:8,2:0.2:10\n'
    return ('##fileformat=VCFv4.2\n'
            '##source=Mutect2\n'
            '##contig=<ID=chrSynthetic,length=8>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
            + '\t'.join(samples) + '\n' + body)


def parse(tmp_path, text, **kwargs):
    from ontseq_platform.gatk_adapter.vcf import extract_candidates
    path = tmp_path / 'input.vcf'
    path.write_text(text)
    return extract_candidates(path, run_id='synthetic-run-001',
                              sample_id='SYNTHETIC_T', **kwargs)


def test_disabled_is_not_run_even_when_files_missing(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    Path(config_dict['tumor']['bam']['path']).unlink()
    result = run_mutect2(config(config_dict), tmp_path / 'out')
    assert result.status == 'NOT_RUN'
    assert not (tmp_path / 'out').exists()
    assert not result.clinically_reportable


def test_enabled_without_ont_opt_in_is_not_run(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    config_dict['enabled'] = True
    result = run_mutect2(config(config_dict), tmp_path / 'out')
    assert result.status == 'NOT_RUN'
    assert 'experimental' in ' '.join(result.messages).lower()


@pytest.mark.parametrize('field', ['surprise', 'bqsr', 'mark_duplicates', 'orientation_bias'])
def test_unknown_options_rejected(config_dict, field):
    from pydantic import ValidationError
    config_dict[field] = True
    with pytest.raises(ValidationError):
        config(config_dict)


def test_paired_requires_normal(config_dict):
    from pydantic import ValidationError
    config_dict['mode'] = 'tumor_normal'
    with pytest.raises(ValidationError):
        config(config_dict)


def test_tumor_only_rejects_normal(config_dict, tmp_path):
    from pydantic import ValidationError
    with_normal(config_dict, tmp_path)['mode'] = 'tumor_only'
    with pytest.raises(ValidationError):
        config(config_dict)


def test_same_sample_id_rejected(config_dict, tmp_path):
    from pydantic import ValidationError
    with_normal(config_dict, tmp_path)['normal']['sample_id'] = 'SYNTHETIC_T'
    with pytest.raises(ValidationError):
        config(config_dict)


def test_same_bam_digest_rejected(config_dict, tmp_path):
    from pydantic import ValidationError
    with_normal(config_dict, tmp_path)['normal']['bam']['sha256'] = config_dict['tumor']['bam']['sha256']
    with pytest.raises(ValidationError):
        config(config_dict)


def test_reference_mismatch_rejected(config_dict):
    from pydantic import ValidationError
    config_dict['tumor']['reference_sha256'] = 'f' * 64
    with pytest.raises(ValidationError):
        config(config_dict)


def test_population_resource_reference_mismatch_rejected(config_dict):
    from pydantic import ValidationError
    config_dict['germline_resource']['reference_sha256'] = 'f' * 64
    with pytest.raises(ValidationError):
        config(config_dict)


def test_contamination_requires_population_sites(config_dict):
    from pydantic import ValidationError
    config_dict['estimate_contamination'] = True
    with pytest.raises(ValidationError):
        config(config_dict)


def test_wrong_version_lock_rejected(config_dict):
    from pydantic import ValidationError
    config_dict['runtime']['expected_version'] = 'latest'
    with pytest.raises(ValidationError):
        config(config_dict)


@pytest.mark.parametrize('value', ['https://example.org/a.bam', 'relative.bam', '/tmp/a\nb.bam'])
def test_nonlocal_or_nonabsolute_paths_rejected(config_dict, value):
    from pydantic import ValidationError
    config_dict['tumor']['bam']['path'] = value
    with pytest.raises(ValidationError):
        config(config_dict)


def test_checksum_drift_fails_preflight(config_dict):
    from ontseq_platform.gatk_adapter.core import preflight
    Path(config_dict['tumor']['bam']['path']).write_text('altered\n')
    with pytest.raises(ValueError, match='checksum'):
        preflight(config(config_dict))


def test_sidecar_path_must_match(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import preflight
    config_dict['tumor']['bai'] = locked(tmp_path / 'other.bai')
    with pytest.raises(ValueError, match='index'):
        preflight(config(config_dict))


def test_reference_dictionary_length_mismatch(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import preflight
    config_dict['reference']['dictionary'] = locked(
        tmp_path / 'reference.dict', '@SQ\tSN:chrSynthetic\tLN:9\n')
    with pytest.raises(ValueError, match='dictionary'):
        preflight(config(config_dict))


def test_pon_must_be_assay_matched(config_dict, tmp_path):
    from pydantic import ValidationError
    resource = dict(config_dict['germline_resource'])
    resource.update(kind='ont_pon', assay_id='different-assay')
    config_dict['panel_of_normals'] = resource
    with pytest.raises(ValidationError):
        config(config_dict)


def test_tumor_only_plan_never_invents_a_normal(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import build_plan
    plan = build_plan(config(config_dict), tmp_path / 'out')
    calling = next(s for s in plan if s.name == 'mutect2')
    assert calling.argv.count('-I') == 1
    assert '-normal' not in calling.argv
    assert '--germline-resource' in calling.argv
    filtering = next(s for s in plan if s.name == 'filter')
    assert '--stats' in filtering.argv
    assert '--f1r2-tar-gz' not in calling.argv
    assert not any('BaseRecalibrator' in s.argv or 'MarkDuplicates' in s.argv for s in plan)


def test_tumor_normal_plan_names_normal_from_manifest(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import build_plan
    plan = build_plan(config(with_normal(config_dict, tmp_path)), tmp_path / 'out')
    calling = next(s for s in plan if s.name == 'mutect2')
    assert calling.argv.count('-I') == 2
    assert calling.argv[calling.argv.index('-normal') + 1] == 'SYNTHETIC_N'


def test_contamination_plan_uses_intersection_and_normal(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import build_plan
    with_normal(config_dict, tmp_path)
    config_dict['estimate_contamination'] = True
    config_dict['contamination_sites'] = config_dict['germline_resource']
    config_dict['intervals'] = locked(tmp_path / 'targets.bed', 'chrSynthetic\t0\t8\n')
    plan = build_plan(config(config_dict), tmp_path / 'out')
    pileups = [s for s in plan if 'pileup' in s.name]
    assert len(pileups) == 2
    assert all('INTERSECTION' in s.argv for s in pileups)
    estimation = next(s for s in plan if s.name == 'contamination')
    assert '--matched-normal' in estimation.argv
    filtering = next(s for s in plan if s.name == 'filter')
    assert '--contamination-table' in filtering.argv
    assert '--tumor-segmentation' in filtering.argv


def test_paths_with_spaces_are_single_arguments(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import build_plan
    config_dict['intervals'] = locked(tmp_path / 'targets with spaces.bed')
    plan = build_plan(config(config_dict), tmp_path / 'out with spaces')
    calling = next(s for s in plan if s.name == 'mutect2')
    assert str(tmp_path / 'targets with spaces.bed') in calling.argv


def test_vcf_coordinate_conversion_and_af(tmp_path):
    records = parse(tmp_path, vcf_text())
    record = records[0]
    assert (record.start0, record.end0, record.native_pos1) == (1, 2, 2)
    assert (record.ref, record.alt, record.variant_type) == ('C', 'T', 'SNV')
    assert record.tumor_af == 0.2
    assert (record.ref_depth, record.alt_depth, record.depth) == (8, 2, 10)
    assert record.native_site_filter == 'PASS'
    assert not record.clinically_reportable
    assert not record.left_aligned


def test_missing_af_not_replaced_with_ad_fraction(tmp_path):
    text = vcf_text().replace('GT:AD:AF:DP\t0/1:8,2:0.2:10', 'GT:AD:DP\t0/1:8,2:10')
    assert parse(tmp_path, text)[0].tumor_af is None


def test_filters_preserved_including_not_applied(tmp_path):
    for value in ['.', 'germline', 'weak_evidence;strand_bias']:
        text = vcf_text().replace('\tPASS\t', f'\t{value}\t')
        assert parse(tmp_path, text)[0].native_site_filter == value


def test_multiallelic_fields_are_allele_specific(tmp_path):
    body = ('chrSynthetic\t3\t.\tG\tA,T\t.\tPASS\tAS_FilterStatus=SITE|weak_evidence'
            '\tGT:AD:AF:DP\t1/2:70,20,10:0.2,0.1:100\n')
    records = parse(tmp_path, vcf_text(body=body))
    assert [r.alt for r in records] == ['A', 'T']
    assert [r.tumor_af for r in records] == [0.2, 0.1]
    assert [r.alt_depth for r in records] == [20, 10]
    assert [r.allele_index for r in records] == [1, 2]
    assert len({r.evidence_id for r in records}) == 2
    assert all(r.native_allele_filter_status == 'SITE|weak_evidence' for r in records)


@pytest.mark.parametrize('af', ['NaN', 'Inf', '-0.1', '1.1'])
def test_invalid_af_rejected(tmp_path, af):
    with pytest.raises(ValueError):
        parse(tmp_path, vcf_text().replace('8,2:0.2:10', f'8,2:{af}:10'))


@pytest.mark.parametrize('suffix', ['0/1:8:0.2:10', '0/1:8,2:0.2,0.1:10', '0/1:-8,2:0.2:10'])
def test_bad_format_cardinality_and_negative_counts_rejected(tmp_path, suffix):
    with pytest.raises(ValueError):
        parse(tmp_path, vcf_text().replace('0/1:8,2:0.2:10', suffix))


def test_wrong_vcf_sample_rejected(tmp_path):
    with pytest.raises(ValueError, match='sample'):
        parse(tmp_path, vcf_text(samples=('WRONG_SAMPLE',)))


def test_extra_vcf_sample_rejected(tmp_path):
    with pytest.raises(ValueError, match='sample'):
        parse(tmp_path, vcf_text(samples=('SYNTHETIC_T', 'UNEXPECTED')))


def test_paired_vcf_sample_order_not_assumed(tmp_path):
    body = ('chrSynthetic\t2\t.\tC\tT\t.\tPASS\t.\tGT:AD:AF:DP'
            '\t0/0:19,1:0.05:20\t0/1:8,2:0.2:10\n')
    records = parse(tmp_path, vcf_text(('SYNTHETIC_N', 'SYNTHETIC_T'), body),
                    normal_sample_id='SYNTHETIC_N')
    assert records[0].tumor_af == 0.2
    assert records[0].normal_af == 0.05


def test_indel_native_representation_is_preserved(tmp_path):
    body = 'chrSynthetic\t2\t.\tCG\tC\t.\tPASS\t.\tGT:AD:AF:DP\t0/1:8,2:0.2:10\n'
    record = parse(tmp_path, vcf_text(body=body))[0]
    assert (record.start0, record.end0, record.ref, record.alt) == (1, 3, 'CG', 'C')
    assert record.variant_type == 'DEL'


def test_header_only_is_not_a_biological_negative(tmp_path):
    assert parse(tmp_path, vcf_text(body='')) == ()


def test_missing_header_rejected(tmp_path):
    with pytest.raises(ValueError, match='header'):
        parse(tmp_path, 'not a VCF\n')


def test_symbolic_allele_rejected_without_silent_drop(tmp_path):
    with pytest.raises(ValueError, match='allele'):
        parse(tmp_path, vcf_text().replace('\tC\tT\t', '\tC\t<DEL>\t'))


def test_gzipped_vcf_supported(tmp_path):
    from ontseq_platform.gatk_adapter.vcf import extract_candidates
    path = tmp_path / 'input.vcf.gz'
    with gzip.open(path, 'wt') as handle:
        handle.write(vcf_text())
    assert len(extract_candidates(path, run_id='synthetic-run-001',
                                  sample_id='SYNTHETIC_T')) == 1


class SyntheticExecutor:
    """Explicit external-process test double, NEVER biological/tool evidence."""
    def __init__(self, *, fail_at=None, wrong_sample=False, wrong_version=False,
                 empty_vcf=False, omit_stats=False):
        self.calls = []
        self.fail_at = fail_at
        self.wrong_sample = wrong_sample
        self.wrong_version = wrong_version
        self.empty_vcf = empty_vcf
        self.omit_stats = omit_stats

    def __call__(self, step, workdir, timeout):
        self.calls.append(step.name)
        out = workdir / f'{step.name}.stdout.log'
        err = workdir / f'{step.name}.stderr.log'
        out.write_text('synthetic process output\n')
        err.write_text('')
        if self.fail_at == step.name:
            return 9
        if step.name == 'version':
            out.write_text('The Genome Analysis Toolkit (GATK) v'
                           + ('4.7.0.0' if self.wrong_version else '4.6.2.0'))
        elif step.name == 'java_version':
            err.write_text('openjdk version "17.0.1"\n')
        elif step.name.endswith('_sample'):
            name = 'SYNTHETIC_N' if step.name.startswith('normal') else 'SYNTHETIC_T'
            step.outputs[0].write_text('WRONG' if self.wrong_sample else name)
        else:
            for output in step.outputs:
                if self.omit_stats and output.name.endswith('.vcf.gz.stats'):
                    continue
                if output.name.endswith('.vcf.gz'):
                    paired = '-normal' in step.argv or 'normal_sample' in self.calls
                    text = vcf_text(body='' if self.empty_vcf else None)
                    if paired:
                        body = ('chrSynthetic\t2\t.\tC\tT\t.\tPASS\t.\tGT:AD:AF:DP'
                                '\t0/1:8,2:0.2:10\t0/0:20,0:0.0:20\n')
                        text = vcf_text(('SYNTHETIC_T', 'SYNTHETIC_N'), body)
                    with gzip.open(output, 'wt') as handle:
                        handle.write(text)
                elif output.name == 'contamination.table':
                    output.write_text('sample\tcontamination\terror\nSYNTHETIC_T\t0.01\t0.001\n')
                else:
                    output.write_text('synthetic output\n')
        return 0


def enabled(data):
    data['enabled'] = True
    data['allow_experimental_ont'] = True
    return config(data)


def test_success_is_candidate_only_with_provenance(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    executor = SyntheticExecutor()
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=executor)
    assert result.status == 'COMPLETED'
    assert len(result.candidates) == 1
    assert not result.clinically_reportable
    assert not result.analytically_validated
    assert result.no_biological_negative_claim
    assert result.execution_evidence == 'TEST_DOUBLE'
    assert result.input_sha256s['tumor.bam'] == config_dict['tumor']['bam']['sha256']
    assert 'filtered.vcf.gz' in result.output_sha256s
    assert json.loads((tmp_path / 'out' / 'result.json').read_text())['status'] == 'COMPLETED'


@pytest.mark.parametrize('failure', ['version', 'tumor_validate', 'tumor_sample', 'mutect2', 'filter'])
def test_failed_step_stops_and_never_becomes_negative(config_dict, tmp_path, failure):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    executor = SyntheticExecutor(fail_at=failure)
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=executor)
    assert result.status == 'FAILED'
    assert not result.candidates
    assert executor.calls[-1] == failure
    assert (tmp_path / 'out' / 'result.json').exists()


def test_wrong_binary_version_blocks_calling(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    executor = SyntheticExecutor(wrong_version=True)
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=executor)
    assert result.status == 'FAILED'
    assert 'mutect2' not in executor.calls


def test_wrong_bam_sample_blocks_calling(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    executor = SyntheticExecutor(wrong_sample=True)
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=executor)
    assert result.status == 'FAILED'
    assert 'mutect2' not in executor.calls


def test_missing_mutect_stats_blocks_filtering(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    executor = SyntheticExecutor(omit_stats=True)
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=executor)
    assert result.status == 'FAILED'
    assert 'filter' not in executor.calls


def test_empty_valid_vcf_not_negative(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    result = run_mutect2(enabled(config_dict), tmp_path / 'out',
                         executor=SyntheticExecutor(empty_vcf=True))
    assert result.status == 'COMPLETED'
    assert result.candidates == ()
    assert result.no_biological_negative_claim


def test_output_directory_not_reused(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    out = tmp_path / 'out'
    out.mkdir()
    marker = out / 'old-result.txt'
    marker.write_text('must remain untouched')
    with pytest.raises(FileExistsError):
        run_mutect2(enabled(config_dict), out, executor=SyntheticExecutor())
    assert marker.read_text() == 'must remain untouched'


def test_input_drift_stops_before_native_execution(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    cfg = enabled(config_dict)
    Path(config_dict['tumor']['bam']['path']).write_text('changed')
    executor = SyntheticExecutor()
    result = run_mutect2(cfg, tmp_path / 'out', executor=executor)
    assert result.status == 'FAILED'
    assert executor.calls == []


def test_paired_and_contamination_execution(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    with_normal(config_dict, tmp_path)
    config_dict['estimate_contamination'] = True
    config_dict['contamination_sites'] = config_dict['germline_resource']
    executor = SyntheticExecutor()
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=executor)
    assert result.status == 'COMPLETED'
    assert 'normal_pileup' in executor.calls
    assert 'contamination' in executor.calls


def test_timeout_recorded_as_failed(config_dict, tmp_path):
    import subprocess
    from ontseq_platform.gatk_adapter.core import run_mutect2
    def timed_out(step, workdir, timeout):
        raise subprocess.TimeoutExpired(list(step.argv), timeout)
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=timed_out)
    assert result.status == 'FAILED'
    assert result.candidates == ()


def test_cli_plan_never_executes_and_reports_not_run(config_dict, tmp_path):
    import os
    import subprocess
    import sys
    config_path = tmp_path / 'configuration.json'
    config_path.write_text(json.dumps(config_dict))
    command = [sys.executable, '-m', 'ontseq_platform.gatk_adapter', 'plan',
               '--config', str(config_path), '--output-dir', str(tmp_path / 'never-created')]
    completed = subprocess.run(command, capture_output=True, text=True, env=os.environ.copy())
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)['status'] == 'NOT_RUN'
    assert not (tmp_path / 'never-created').exists()


def test_cli_run_requires_enablement(config_dict, tmp_path):
    import os
    import subprocess
    import sys
    config_path = tmp_path / 'configuration.json'
    config_path.write_text(json.dumps(config_dict))
    completed = subprocess.run([
        sys.executable, '-m', 'ontseq_platform.gatk_adapter', 'run',
        '--config', str(config_path), '--output-dir', str(tmp_path / 'never-created'),
    ], capture_output=True, text=True, env=os.environ.copy())
    assert completed.returncode == 2
    assert json.loads(completed.stdout)['status'] == 'NOT_RUN'
    assert not (tmp_path / 'never-created').exists()


def test_manifest_builder_hashes_actual_local_inputs(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.cli import lock_config
    built = lock_config(
        run_id='synthetic-run-002', assay_id='synthetic-assay',
        reference_id='synthetic-ref', genome_build='synthetic',
        reference=Path(config_dict['reference']['fasta']['path']),
        tumor_bam=Path(config_dict['tumor']['bam']['path']),
        tumor_sample='SYNTHETIC_T',
        germline=Path(config_dict['germline_resource']['vcf']['path']),
        gatk_jar=Path(config_dict['runtime']['gatk_jar']['path']),
    )
    assert built.tumor.bam.sha256 == config_dict['tumor']['bam']['sha256']
    assert not built.enabled and not built.allow_experimental_ont


def test_contamination_nan_fails_closed(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    config_dict['estimate_contamination'] = True
    config_dict['contamination_sites'] = config_dict['germline_resource']
    fake = SyntheticExecutor()
    def execute(step, workdir, timeout):
        code = fake(step, workdir, timeout)
        if step.name == 'contamination':
            (workdir / 'contamination.table').write_text(
                'sample\tcontamination\terror\nSYNTHETIC_T\tNaN\t0.001\n')
        return code
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=execute)
    assert result.status == 'FAILED'
    assert 'filter' not in fake.calls


def test_input_mutation_during_run_rejects_candidates(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    fake = SyntheticExecutor()
    def execute(step, workdir, timeout):
        code = fake(step, workdir, timeout)
        if step.name == 'filter':
            Path(config_dict['tumor']['bam']['path']).write_text('changed during run')
        return code
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=execute)
    assert result.status == 'FAILED'
    assert not result.candidates


def test_unregistered_secondary_index_rejected(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import preflight
    (tmp_path / 'tumor.bai').write_text('unregistered-second-index')
    with pytest.raises(ValueError, match='index'):
        preflight(config(config_dict))


def test_mutect_result_cannot_be_clinically_released(config_dict):
    from pydantic import ValidationError
    from ontseq_platform.gatk_adapter.contracts import Mutect2Result
    from ontseq_platform.gatk_adapter.core import config_sha256
    with pytest.raises(ValidationError):
        Mutect2Result(
            run_id='synthetic-run-001', sample_id='SYNTHETIC_T', mode='tumor_only',
            status='COMPLETED', clinically_reportable=True,
            config_sha256=config_sha256(config(config_dict)),
        )


def test_locked_file_and_config_json_roundtrip(config_dict):
    from ontseq_platform.gatk_adapter.contracts import LockedFile, Mutect2Config
    original = config(config_dict)
    encoded = original.model_dump_json()
    assert Mutect2Config.model_validate_json(encoded) == original
    artifact = original.reference.fasta
    assert LockedFile.model_validate_json(artifact.model_dump_json()) == artifact


def test_preflight_failure_does_not_claim_any_process_execution(config_dict, tmp_path):
    from ontseq_platform.gatk_adapter.core import run_mutect2
    fake = SyntheticExecutor()
    Path(config_dict['tumor']['bam']['path']).write_text('changed before start')
    result = run_mutect2(enabled(config_dict), tmp_path / 'out', executor=fake)
    assert result.status == 'FAILED'
    assert fake.calls == []
    assert result.execution_evidence == 'NONE'


def test_native_process_transport_is_shell_free_and_scrubs_java_environment(tmp_path, monkeypatch):
    import os
    import sys
    from ontseq_platform.gatk_adapter.core import Step, _native_execute
    marker = tmp_path / 'must-not-exist'
    payload = f'; touch {marker}'
    for key in ('JAVA_TOOL_OPTIONS', 'JDK_JAVA_OPTIONS', '_JAVA_OPTIONS', 'CLASSPATH'):
        monkeypatch.setenv(key, 'must-not-propagate')
    program = (
        'import json, os, sys; '
        'print(json.dumps({"argument": sys.argv[1], "env": '
        '{k: os.environ.get(k) for k in '
        '["JAVA_TOOL_OPTIONS", "JDK_JAVA_OPTIONS", "_JAVA_OPTIONS", "CLASSPATH"]}}))'
    )
    # This is a harmless Python subprocess transport test, NOT a native GATK test.
    result = _native_execute(Step('transport', (sys.executable, '-c', program, payload)),
                             tmp_path, 10)
    assert result == 0
    output = json.loads((tmp_path / 'transport.stdout.log').read_text())
    assert output['argument'] == payload
    assert set(output['env'].values()) == {None}
    assert not marker.exists()
    if os.name == 'posix':
        assert (tmp_path / 'transport.stdout.log').stat().st_mode & 0o777 == 0o600


def test_cli_lock_writes_disabled_json_that_can_be_read(config_dict, tmp_path):
    import os
    import subprocess
    import sys
    from ontseq_platform.gatk_adapter.contracts import Mutect2Config
    output = tmp_path / 'locked-cli.json'
    command = [sys.executable, '-m', 'ontseq_platform.gatk_adapter', 'lock',
               '--run-id', 'synthetic-cli-001', '--assay-id', 'synthetic-assay',
               '--reference-id', 'synthetic-ref', '--genome-build', 'synthetic',
               '--reference', config_dict['reference']['fasta']['path'],
               '--tumor-bam', config_dict['tumor']['bam']['path'],
               '--tumor-sample', 'SYNTHETIC_T',
               '--germline', config_dict['germline_resource']['vcf']['path'],
               '--gatk-jar', config_dict['runtime']['gatk_jar']['path'],
               '--output', str(output)]
    completed = subprocess.run(command, capture_output=True, text=True, env=os.environ.copy())
    assert completed.returncode == 0, completed.stderr
    restored = Mutect2Config.model_validate_json(output.read_text())
    assert not restored.enabled and not restored.allow_experimental_ont
    assert restored.tumor.bam.sha256 == config_dict['tumor']['bam']['sha256']


def test_isolated_schemas_match_current_contracts():
    from ontseq_platform.gatk_adapter.contracts import (
        Mutect2Config, Mutect2Result, SmallVariantCandidate,
    )
    directory = Path(__file__).resolve().parents[1] / 'schemas' / 'gatk_adapter'
    for name, model in [('config', Mutect2Config), ('result', Mutect2Result),
                        ('candidate', SmallVariantCandidate)]:
        schema = model.model_json_schema()
        schema['$schema'] = 'https://json-schema.org/draft/2020-12/schema'
        assert json.loads((directory / (name + '.schema.json')).read_text()) == schema
