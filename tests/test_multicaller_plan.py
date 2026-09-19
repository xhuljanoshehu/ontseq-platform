from __future__ import annotations

import unittest


class MultiCallerCatalogTests(unittest.TestCase):
    def test_catalog_has_eight_distinct_callers(self) -> None:
        from ontseq_platform.multicaller_catalog import CALLER_CATALOG

        self.assertEqual(
            set(CALLER_CATALOG),
            {
                "qdnaseq_ace",
                "ont_spectre",
                "ichorcna",
                "sniffles2",
                "cutesv",
                "savana",
                "severus",
                "wakhan",
            },
        )


class MultiCallerCatalogContractTests(unittest.TestCase):
    def test_catalog_separates_existing_adapters_from_planned_adapters(self) -> None:
        from ontseq_platform.multicaller_catalog import CALLER_CATALOG

        self.assertEqual(
            {key for key, entry in CALLER_CATALOG.items() if entry.existing_adapter},
            {"qdnaseq_ace", "sniffles2", "cutesv"},
        )

    def test_breakpoint_callers_are_not_absolute_copy_number_callers(self) -> None:
        from ontseq_platform.multicaller_catalog import CALLER_CATALOG

        for key in ("sniffles2", "cutesv", "severus"):
            self.assertNotIn("cnv", CALLER_CATALOG[key].analyses)
        self.assertNotIn("purity_ploidy", CALLER_CATALOG["ont_spectre"].analyses)
        self.assertIn("allelic_cn", CALLER_CATALOG["wakhan"].analyses)

    def test_catalog_is_immutable(self) -> None:
        from dataclasses import FrozenInstanceError

        from ontseq_platform.multicaller_catalog import CALLER_CATALOG

        with self.assertRaises(TypeError):
            CALLER_CATALOG["invented"] = CALLER_CATALOG["severus"]
        with self.assertRaises(FrozenInstanceError):
            CALLER_CATALOG["severus"].label = "changed"


def _sha(label: str) -> str:
    import hashlib

    return hashlib.sha256(label.encode()).hexdigest()


def _context(**updates):
    payload = {
        "sample_id": "SYNTHETIC_TUMOUR_001",
        "genome_build": "GRCh38",
        "reference_sha256": _sha("reference"),
        "bam_sha256": _sha("bam"),
        "assessability_mask_sha256": _sha("mask"),
        "platform": "ONT",
        "data_basis": "lcwgs_genome_wide",
        "coverage_x": 12.0,
        "coverage_definition": "synthetic mean autosomal depth",
        "artifacts": [],
    }
    payload.update(updates)
    if payload["coverage_x"] is None and "coverage_definition" not in updates:
        payload["coverage_definition"] = None
    return payload


def _artifact(kind: str, *, resource: bool = False, **updates):
    payload = {
        "kind": kind,
        "sha256": _sha(kind),
        "genome_build": "GRCh38",
        "reference_sha256": _sha("reference"),
        "sample_id": None if resource else "SYNTHETIC_TUMOUR_001",
        "source_input_sha256": None if resource else _sha("bam"),
    }
    payload.update(updates)
    return payload


def _selection(caller="qdnaseq_ace", analyses=None, mode="tumour_only", **updates):
    payload = {
        "caller_id": caller,
        "analyses": analyses or ["cnv"],
        "mode": mode,
        "policy": {
            "caller_id": caller,
            "genome_build": "GRCh38",
            "platform": "ONT",
            "data_basis": "lcwgs_genome_wide",
            "allowed_modes": [mode],
            "minimum_coverage_x": 5.0,
            "coverage_definition": "synthetic mean autosomal depth",
            "reference_sha256": _sha("reference"),
            "assessability_mask_sha256": _sha("mask"),
        },
        "runtime": {
            "version": "1.0.0",
            "runtime_sha256": _sha("runtime"),
            "adapter_sha256": _sha("adapter"),
            "parameters_sha256": _sha("parameters"),
        },
    }
    payload.update(updates)
    return payload


def _request(*selections, context=None):
    return {
        "study_registration_sha256": _sha("study"),
        "context": context or _context(),
        "selections": list(selections),
    }


def _plan(payload):
    from ontseq_platform.multicaller_plan import MultiCallerRequest, plan_multicaller

    return plan_multicaller(MultiCallerRequest.model_validate(payload))


def _lane(plan, caller="qdnaseq_ace"):
    return next(item for item in plan.lanes if item.caller_id == caller)


class MultiCallerRoutingTests(unittest.TestCase):
    def test_existing_adapter_is_planned_not_executed(self) -> None:
        plan = _plan(_request(_selection()))
        lane = _lane(plan)
        self.assertEqual(lane.disposition, "PLANNED")
        self.assertTrue(lane.input_eligible)
        self.assertTrue(lane.existing_adapter)
        self.assertEqual(lane.execution_status, "NOT_RUN")
        self.assertFalse(plan.execution_enabled)
        self.assertTrue(plan.research_only)
        self.assertTrue(plan.human_review_required)
        self.assertTrue(plan.retain_all_evidence)
        self.assertEqual(plan.comparison_policy, "NO_VOTING")

    def test_all_eight_callers_remain_visible(self) -> None:
        plan = _plan(_request(_selection()))
        self.assertEqual(len(plan.lanes), 8)
        self.assertEqual(sum(row.disposition == "NOT_SELECTED" for row in plan.lanes), 7)
        self.assertTrue(all(row.execution_status == "NOT_RUN" for row in plan.lanes))

    def test_disabled_caller_is_not_an_input_failure_or_a_negative(self) -> None:
        plan = _plan(_request(_selection(enabled=False, policy=None, runtime=None)))
        self.assertEqual(_lane(plan).disposition, "NOT_SELECTED")
        self.assertEqual(_lane(plan).blockers, ())
        self.assertFalse(_lane(plan).input_eligible)

    def test_unknown_caller_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _plan(_request(_selection("spectre")))

    def test_duplicate_caller_request_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _plan(_request(_selection(), _selection()))

    def test_duplicate_analysis_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _plan(_request(_selection(analyses=["cnv", "cnv"])))

    def test_unsupported_analysis_is_blocked(self) -> None:
        lane = _lane(_plan(_request(_selection("sniffles2"))), "sniffles2")
        self.assertIn("ANALYSIS_NOT_SUPPORTED", lane.blockers)

    def test_missing_policy_and_runtime_are_both_visible(self) -> None:
        lane = _lane(_plan(_request(_selection(policy=None, runtime=None))))
        self.assertIn("STUDY_POLICY_MISSING", lane.blockers)
        self.assertIn("RUNTIME_LOCK_MISSING", lane.blockers)
        self.assertEqual(lane.disposition, "BLOCKED")

    def test_unknown_coverage_is_not_imputed(self) -> None:
        plan = _plan(_request(_selection(), context=_context(coverage_x=None)))
        self.assertIn("COVERAGE_UNKNOWN", _lane(plan).blockers)
        self.assertIsNone(plan.request.context.coverage_x)
        self.assertIsNone(plan.request.context.tumor_fraction)

    def test_nonfinite_negative_and_boolean_coverage_are_rejected(self) -> None:
        for value in (float("nan"), float("inf"), -1.0, True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                _plan(_request(_selection(), context=_context(coverage_x=value)))

    def test_coverage_threshold_is_prospective_and_inclusive(self) -> None:
        below = _plan(_request(_selection(), context=_context(coverage_x=4.99)))
        equal = _plan(_request(_selection(), context=_context(coverage_x=5.0)))
        self.assertIn("COVERAGE_BELOW_REGISTERED_MINIMUM", _lane(below).blockers)
        self.assertTrue(_lane(equal).input_eligible)

    def test_policy_identity_must_match_caller_assay_build_and_platform(self) -> None:
        for field, value in (
            ("caller_id", "savana"),
            ("genome_build", "GRCh37"),
            ("data_basis", "adaptive_sampling_on_target"),
            ("platform", "ILLUMINA"),
        ):
            selection = _selection()
            selection["policy"][field] = value
            with self.subTest(field=field):
                self.assertIn(
                    "STUDY_POLICY_IDENTITY_MISMATCH",
                    _lane(_plan(_request(selection))).blockers,
                )

    def test_mode_outside_registered_policy_is_blocked(self) -> None:
        selection = _selection()
        selection["policy"]["allowed_modes"] = ["single_sample"]
        self.assertIn("MODE_NOT_REGISTERED", _lane(_plan(_request(selection))).blockers)

    def test_moving_runtime_versions_are_rejected(self) -> None:
        for version in ("latest", "main", "master", "1.*", ">=1.0", "", " latest "):
            selection = _selection()
            selection["runtime"]["version"] = version
            with self.subTest(version=version), self.assertRaises(ValueError):
                _plan(_request(selection))

    def test_ref_or_sample_mismatch_is_rejected_before_routing(self) -> None:
        for update in (
            {"reference_sha256": _sha("other-reference")},
            {"genome_build": "GRCh37"},
            {"sample_id": "SYNTHETIC_OTHER"},
            {"source_input_sha256": _sha("other-bam")},
        ):
            with self.subTest(update=update), self.assertRaises(ValueError):
                _plan(_request(context=_context(artifacts=[_artifact("snv_vcf", **update)])))

    def test_duplicate_artifact_kind_is_rejected(self) -> None:
        artifact = _artifact("snv_vcf")
        with self.assertRaises(ValueError):
            _plan(_request(context=_context(artifacts=[artifact, artifact])))

    def test_tumour_cannot_be_its_own_normal(self) -> None:
        for update in (
            {"normal_sample_id": "SYNTHETIC_TUMOUR_001", "normal_bam_sha256": _sha("normal")},
            {"normal_sample_id": "SYNTHETIC_NORMAL", "normal_bam_sha256": _sha("bam")},
        ):
            with self.subTest(update=update), self.assertRaises(ValueError):
                _plan(_request(context=_context(**update)))

    def test_paired_somatic_request_does_not_silently_fall_back(self) -> None:
        selection = _selection("severus", ["sv"], "tumour_normal")
        lane = _lane(_plan(_request(selection)), "severus")
        self.assertIn("MATCHED_NORMAL_MISSING", lane.blockers)
        self.assertEqual(lane.mode, "tumour_normal")

    def test_paired_normal_depth_is_checked_independently(self) -> None:
        selection = _selection("severus", ["sv"], "tumour_normal")
        context = _context(normal_sample_id="SYNTHETIC_NORMAL", normal_bam_sha256=_sha("normal"))
        lane = _lane(_plan(_request(selection, context=context)), "severus")
        self.assertIn("NORMAL_COVERAGE_POLICY_MISSING", lane.blockers)
        self.assertIn("NORMAL_COVERAGE_UNKNOWN", lane.blockers)
        selection["policy"]["minimum_normal_coverage_x"] = 8.0
        selection["policy"]["normal_coverage_definition"] = "synthetic normal mean depth"
        context["normal_coverage_definition"] = "synthetic normal mean depth"
        context["normal_coverage_x"] = 7.0
        lane = _lane(_plan(_request(selection, context=context)), "severus")
        self.assertIn("NORMAL_COVERAGE_BELOW_REGISTERED_MINIMUM", lane.blockers)
        context["normal_coverage_x"] = 8.0
        lane = _lane(_plan(_request(selection, context=context)), "severus")
        self.assertTrue(lane.input_eligible)
        self.assertEqual(lane.disposition, "ADAPTER_PENDING")

    def test_tumour_only_sv_is_not_labeled_confirmed_somatic(self) -> None:
        lane = _lane(_plan(_request(_selection("severus", ["sv"]))), "severus")
        self.assertIn("UNPAIRED_SV_IS_CANDIDATE_EVIDENCE", lane.warnings)
        self.assertEqual(lane.execution_status, "NOT_RUN")

    def test_adaptive_sampling_cn_is_not_silently_treated_as_wgs(self) -> None:
        for basis in ("adaptive_sampling_off_target", "adaptive_sampling_on_target"):
            selection = _selection()
            selection["policy"]["data_basis"] = basis
            lane = _lane(_plan(_request(selection, context=_context(data_basis=basis))))
            self.assertIn("AS_GENOMEWIDE_CN_NOT_QUALIFIED", lane.blockers)

    def test_adaptive_sampling_sv_keeps_a_separate_registered_scope(self) -> None:
        selection = _selection("sniffles2", ["sv"])
        selection["policy"]["data_basis"] = "adaptive_sampling_on_target"
        context = _context(data_basis="adaptive_sampling_on_target")
        lane = _lane(_plan(_request(selection, context=context)), "sniffles2")
        self.assertTrue(lane.input_eligible)
        self.assertIn("AS_BREAKPOINT_ASSESSABILITY_REQUIRES_REVIEW", lane.warnings)

    def test_spectre_requires_sample_snvs_coverage_bins_and_registered_bin_size(self) -> None:
        selection = _selection("ont_spectre")
        lane = _lane(_plan(_request(selection)), "ont_spectre")
        for blocker in ("SAMPLE_SNV_REQUIRED", "COVERAGE_BINS_REQUIRED", "BIN_SIZE_UNREGISTERED"):
            self.assertIn(blocker, lane.blockers)
        selection["policy"]["bin_size_bp"] = 1000
        context = _context(
            artifacts=[_artifact("snv_vcf"), _artifact("coverage_bins", bin_size_bp=1000)]
        )
        lane = _lane(_plan(_request(selection, context=context)), "ont_spectre")
        self.assertTrue(lane.input_eligible)
        self.assertFalse(lane.existing_adapter)
        self.assertEqual(lane.disposition, "ADAPTER_PENDING")
        self.assertIn("TUMOUR_DIPLOID_BASELINE_ASSUMPTION_REQUIRES_VALIDATION", lane.warnings)

    def test_spectre_bin_size_mismatch_is_blocked(self) -> None:
        selection = _selection("ont_spectre")
        selection["policy"]["bin_size_bp"] = 1000
        context = _context(
            artifacts=[_artifact("snv_vcf"), _artifact("coverage_bins", bin_size_bp=500)]
        )
        lane = _lane(_plan(_request(selection, context=context)), "ont_spectre")
        self.assertIn("COVERAGE_BIN_SIZE_MISMATCH", lane.blockers)

    def test_savana_sv_and_cn_have_different_snp_requirements(self) -> None:
        sv = _lane(_plan(_request(_selection("savana", ["sv"]))), "savana")
        cn = _lane(_plan(_request(_selection("savana", ["cnv"]))), "savana")
        self.assertTrue(sv.input_eligible)
        self.assertIn("SAVANA_SNP_EVIDENCE_REQUIRED", cn.blockers)
        context = _context(artifacts=[_artifact("snv_vcf")])
        cn = _lane(_plan(_request(_selection("savana", ["cnv"]), context=context)), "savana")
        self.assertTrue(cn.input_eligible)

    def test_wakhan_requires_phased_inputs_and_explicit_segmentation(self) -> None:
        selection = _selection("wakhan", ["allelic_cn"])
        lane = _lane(_plan(_request(selection)), "wakhan")
        for blocker in (
            "PHASED_SNVS_REQUIRED",
            "HAPLOTAGGED_BAM_REQUIRED",
            "SEGMENTATION_UNDECLARED",
        ):
            self.assertIn(blocker, lane.blockers)
        selection["segmentation"] = "change_point"
        context = _context(artifacts=[_artifact("phased_snv_vcf"), _artifact("haplotagged_bam")])
        lane = _lane(_plan(_request(selection, context=context)), "wakhan")
        self.assertTrue(lane.input_eligible)
        self.assertEqual(lane.disposition, "ADAPTER_PENDING")

    def test_sv_guided_wakhan_retains_its_dependency_not_an_independent_vote(self) -> None:
        selection = _selection("wakhan", ["allelic_cn"], segmentation="sv_guided")
        context = _context(artifacts=[_artifact("phased_snv_vcf"), _artifact("haplotagged_bam")])
        lane = _lane(_plan(_request(selection, context=context)), "wakhan")
        self.assertIn("SV_BREAKPOINTS_REQUIRED", lane.blockers)
        context["artifacts"].append(_artifact("sv_breakpoints", producer="severus"))
        lane = _lane(_plan(_request(selection, context=context)), "wakhan")
        self.assertTrue(lane.input_eligible)
        self.assertIn("evidence:severus", lane.dependencies)

    def test_ichorcna_is_an_isolated_cfdna_research_route(self) -> None:
        selection = _selection("ichorcna")
        lane = _lane(_plan(_request(selection)), "ichorcna")
        self.assertIn("OUTSIDE_INITIAL_CALLER_DATA_SCOPE", lane.blockers)
        selection["policy"]["platform"] = "ILLUMINA"
        selection["policy"]["data_basis"] = "cfdna_ulp_wgs"
        context = _context(
            platform="ILLUMINA",
            data_basis="cfdna_ulp_wgs",
            artifacts=[
                _artifact("read_count_wig"),
                _artifact("gc_track", resource=True),
                _artifact("mappability_track", resource=True),
            ],
        )
        lane = _lane(_plan(_request(selection, context=context)), "ichorcna")
        self.assertTrue(lane.input_eligible)
        self.assertEqual(lane.disposition, "ADAPTER_PENDING")

    def test_nonlongread_input_does_not_enable_longread_callers(self) -> None:
        selection = _selection("sniffles2", ["sv"])
        selection["policy"]["platform"] = "ILLUMINA"
        lane = _lane(_plan(_request(selection, context=_context(platform="ILLUMINA"))), "sniffles2")
        self.assertIn("OUTSIDE_INITIAL_PLATFORM_SCOPE", lane.blockers)

    def test_missing_artifact_source_identity_is_rejected(self) -> None:
        artifact = _artifact("coverage_bins", source_input_sha256=None)
        with self.assertRaises(ValueError):
            _plan(_request(context=_context(artifacts=[artifact])))

    def test_shared_population_resource_has_no_fabricated_sample_id(self) -> None:
        artifact = _artifact("population_snps", resource=True)
        plan = _plan(_request(_selection("savana"), context=_context(artifacts=[artifact])))
        self.assertTrue(_lane(plan, "savana").input_eligible)
        artifact["sample_id"] = "SYNTHETIC_TUMOUR_001"
        with self.assertRaises(ValueError):
            _plan(_request(context=_context(artifacts=[artifact])))


class MultiCallerPlanLockTests(unittest.TestCase):
    def test_request_order_does_not_change_the_plan(self) -> None:
        context = _context(artifacts=[_artifact("phased_snv_vcf"), _artifact("haplotagged_bam")])
        selections = [_selection(), _selection("sniffles2", ["sv"])]
        first = _plan(_request(*selections, context=context))
        context["artifacts"].reverse()
        second = _plan(_request(*reversed(selections), context=context))
        self.assertEqual(first.plan_sha256, second.plan_sha256)
        self.assertEqual(first.model_dump(mode="json"), second.model_dump(mode="json"))

    def test_coverage_policy_and_native_fingerprint_changes_change_lock(self) -> None:
        baseline = _plan(_request(_selection()))
        context = _context(coverage_x=13.0)
        self.assertNotEqual(
            baseline.plan_sha256, _plan(_request(_selection(), context=context)).plan_sha256
        )
        selection = _selection()
        selection["policy"]["minimum_coverage_x"] = 6.0
        self.assertNotEqual(baseline.plan_sha256, _plan(_request(selection)).plan_sha256)
        context = _context(bam_sha256=_sha("different-input"))
        self.assertNotEqual(
            baseline.plan_sha256, _plan(_request(_selection(), context=context)).plan_sha256
        )

    def test_locked_plan_rejects_a_tampered_lane(self) -> None:
        from ontseq_platform.multicaller_plan import MultiCallerPlan

        payload = _plan(_request(_selection())).model_dump(mode="json")
        payload["lanes"][0]["input_eligible"] = not payload["lanes"][0]["input_eligible"]
        with self.assertRaises(ValueError):
            MultiCallerPlan.model_validate(payload)

    def test_rehashing_cannot_forge_lane_readiness(self) -> None:
        import hashlib
        import json

        from ontseq_platform.multicaller_plan import MultiCallerPlan

        payload = _plan(_request(_selection())).model_dump(mode="json")
        lane = next(row for row in payload["lanes"] if row["caller_id"] == "wakhan")
        lane["disposition"] = "PLANNED"
        lane["input_eligible"] = True
        payload.pop("plan_sha256")
        payload["plan_sha256"] = hashlib.sha256(
            json.dumps(
                payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
            ).encode()
        ).hexdigest()
        with self.assertRaises(ValueError):
            MultiCallerPlan.model_validate(payload)

    def test_plans_cannot_claim_execution_or_a_winner(self) -> None:
        from ontseq_platform.multicaller_plan import MultiCallerPlan

        for field, value in (
            ("winner", "severus"),
            ("execution_enabled", True),
            ("timestamp", "now"),
        ):
            payload = _plan(_request(_selection())).model_dump(mode="json")
            payload[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                MultiCallerPlan.model_validate(payload)
        payload = _plan(_request(_selection())).model_dump(mode="json")
        payload["lanes"][0]["execution_status"] = "COMPLETED"
        with self.assertRaises(ValueError):
            MultiCallerPlan.model_validate(payload)

    def test_unchecked_model_copy_is_revalidated_at_entry(self) -> None:
        from ontseq_platform.multicaller_plan import MultiCallerRequest, plan_multicaller

        request = MultiCallerRequest.model_validate(_request(_selection()))
        context = request.context.model_copy(update={"coverage_x": float("nan")})
        tampered = request.model_copy(update={"context": context})
        with self.assertRaises(ValueError):
            plan_multicaller(tampered)


class MultiCallerScopeBindingTests(unittest.TestCase):
    def test_reference_mask_and_depth_definition_match_registered_policy(self) -> None:
        for field, value, code in (
            ("reference_sha256", _sha("other-ref"), "STUDY_POLICY_IDENTITY_MISMATCH"),
            ("assessability_mask_sha256", _sha("other-mask"), "STUDY_POLICY_IDENTITY_MISMATCH"),
            ("coverage_definition", "target-only median", "COVERAGE_DEFINITION_MISMATCH"),
        ):
            selection = _selection()
            selection["policy"][field] = value
            with self.subTest(field=field):
                self.assertIn(code, _lane(_plan(_request(selection))).blockers)

    def test_measured_depth_requires_its_definition(self) -> None:
        with self.assertRaises(ValueError):
            _plan(_request(context=_context(coverage_definition=None)))

    def test_unknown_depth_cannot_have_a_measurement_definition(self) -> None:
        with self.assertRaises(ValueError):
            _plan(_request(context=_context(coverage_x=None, coverage_definition="not measured")))

    def test_tumour_fraction_requires_orthogonal_measurement_metadata(self) -> None:
        with self.assertRaises(ValueError):
            _plan(_request(context=_context(tumor_fraction=0.25)))
        context = _context(
            tumor_fraction=0.25,
            tumor_fraction_method="synthetic orthogonal method",
            tumor_fraction_timepoint="same synthetic aliquot",
        )
        plan = _plan(_request(context=context))
        self.assertEqual(plan.request.context.tumor_fraction, 0.25)
        self.assertEqual(plan.request.context.tumor_fraction_timepoint, "same synthetic aliquot")

    def test_normal_snvs_are_bound_to_the_normal_not_to_the_tumour_bam(self) -> None:
        context = _context(
            normal_sample_id="SYNTHETIC_NORMAL",
            normal_bam_sha256=_sha("normal"),
            artifacts=[
                _artifact(
                    "snv_vcf",
                    sample_id="SYNTHETIC_NORMAL",
                    source_input_sha256=_sha("normal"),
                )
            ],
        )
        lane = _lane(_plan(_request(_selection("savana"), context=context)), "savana")
        self.assertTrue(lane.input_eligible)
        context["artifacts"][0]["source_input_sha256"] = _sha("bam")
        with self.assertRaises(ValueError):
            _plan(_request(context=context))


class MultiCallerSchemaTests(unittest.TestCase):
    def test_versioned_schemas_are_current(self) -> None:
        import json
        from pathlib import Path

        from ontseq_platform.multicaller_plan import MultiCallerPlan, MultiCallerRequest

        root = Path(__file__).resolve().parents[1]
        for name, model in (
            ("multicaller-request", MultiCallerRequest),
            ("multicaller-plan", MultiCallerPlan),
        ):
            path = root / "schemas" / f"{name}.schema.json"
            self.assertTrue(path.is_file(), f"Missing public schema: {path.name}")
            self.assertEqual(
                path.read_text(encoding="utf-8"),
                json.dumps(model.model_json_schema(), indent=2, sort_keys=True) + "\n",
            )


if __name__ == "__main__":
    unittest.main()
