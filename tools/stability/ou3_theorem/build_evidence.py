#!/usr/bin/env python3
"""Validate committed theorem status and shipping-source provenance."""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from theorem_status import status_report

REPO=Path(__file__).resolve().parents[3]
PROVENANCE=REPO/"reports/results/ou3_stability/provenance.json"
STATUS=REPO/"reports/results/ou3_stability/theorem-status.json"
RICCATI_STATUS=REPO/"reports/results/ou3_stability/interval-riccati-status.json"

def git_blob_sha(path: Path) -> str:
    data=path.read_bytes(); h=hashlib.sha1(); h.update(f"blob {len(data)}\0".encode("ascii")); h.update(data); return h.hexdigest()

def validate() -> dict:
    provenance=json.loads(PROVENANCE.read_text(encoding="utf-8"))
    committed=json.loads(STATUS.read_text(encoding="utf-8"))
    failures=[]
    sources = provenance["authoritative_shipping_sources"] + provenance.get("operation_lemma_sources", [])
    for row in sources:
        path=REPO/row["path"]
        if not path.is_file(): failures.append(f"missing bound source: {row['path']}"); continue
        actual=git_blob_sha(path)
        if actual!=row["git_blob_sha"]: failures.append(f"source provenance changed: {row['path']} {row['git_blob_sha']} -> {actual}")
    expected=status_report()
    if committed!=expected: failures.append("committed theorem-status.json differs from theorem_status.status_report()")
    riccati=json.loads(RICCATI_STATUS.read_text(encoding="utf-8"))
    if riccati.get("qualification")!="OU3_A21_INTERVAL_RICCATI_V1":
        failures.append("interval Riccati status has wrong qualification")
    if riccati.get("certificate_complete") is True:
        from interval_riccati import load_certificate
        verified=load_certificate(RICCATI_STATUS).verify()
        if not verified["certificate_complete"]:
            failures.append("interval Riccati status claims completion without verified inclusion")
    if str(REPO) not in sys.path:
        sys.path.insert(0,str(REPO))
    from tools.stability.ou3_theorem.lin_matrix_certificate import certificate as matrix_certificate
    from tools.stability.ou3_theorem.word_energy import restricted_service_counterexample
    from tools.stability.ou3_theorem.block_factor_numeric import certificate as block_certificate
    from tools.stability.ou3_theorem.root_covariance_certificate import certificate as root_certificate
    from tools.stability.ou3_theorem.factor_certificates import audit_certificate
    from tools.stability.ou3_theorem.nuisance_upper_certificate import certificate as nuisance_certificate
    from tools.stability.ou3_theorem.sampled_capture_obstruction import witness_certificate
    from tools.stability.ou3_theorem.sampling_fidelity import certificate as sampling_certificate
    from tools.stability.ou3_theorem.corrected_word import certificate as corrected_certificate
    from tools.stability.ou3_theorem.ag_readout import certificate as readout_certificate
    from tools.stability.ou3_theorem.signed_temporal import certificate as signed_certificate
    from tools.stability.ou3_theorem.gyro_bias_projection import certificate as gyro_certificate
    from tools.stability.ou3_theorem.regimes import certificate as regime_certificate
    from tools.stability.ou3_theorem.moving_pivots import certificate as moving_certificate
    from tools.stability.ou3_theorem.stationary_covariance import certificate as stationary_certificate
    from tools.stability.ou3_theorem.regime_continuation_diagnostic import diagnostic as regime_diagnostic
    from tools.stability.ou3_theorem.world_frame import certificate as world_certificate
    from tools.stability.ou3_theorem.aw_covariance_ceiling import certificate as aw_ceiling_certificate
    from tools.stability.ou3_theorem.aw_tracking import certificate as aw_tracking_certificate
    from tools.stability.ou3_theorem.signed_injection import certificate as injection_certificate
    from tools.stability.ou3_theorem.aggregate_floor import certificate as aggregate_certificate
    from tools.stability.ou3_theorem.word_diameter import certificate as diameter_certificate
    from tools.stability.ou3_theorem.imu_two_timescale_certificate import certificate as imu_certificate
    from tools.stability.ou3_theorem.moving_quiet_compatibility import certificate as compatibility_certificate
    from tools.stability.ou3_theorem.planar_service_cell import certificate as planar_cell_certificate
    from tools.stability.ou3_theorem.planar_service_guard import certificate as planar_guard_certificate
    from tools.stability.ou3_theorem.planar_service_frontend_binding import certificate as planar_frontend_binding_certificate
    from tools.stability.ou3_theorem.planar_service_mahony_tube import certificate as planar_mahony_tube_certificate
    from tools.stability.ou3_theorem.planar_generated_port_structure import certificate as planar_generated_port_certificate
    from tools.stability.ou3_theorem.planar_joint_phase_cell import calculate as planar_joint_threshold
    def planar_joint_threshold_certificate():
        return planar_joint_threshold(
            json.loads((STATUS.parent/"planar-service-frechet-diagnostic.json").read_text()),
            planar_mahony_tube_certificate(),
            7.024764605642485)
    from tools.stability.ou3_theorem.mahony_raw_normalization import certificate as raw_norm_certificate
    from tools.stability.ou3_theorem.planar_linked_riccati_mean import certificate as linked_calculus_certificate
    from tools.stability.ou3_theorem.planar_causal_calculus import certificate as causal_calculus_certificate
    from tools.stability.ou3_theorem.planar_mean_covariance_ports import certificate as mean_ports_certificate
    from tools.stability.ou3_theorem.planar_frontend_domain import certificate as frontend_domain_certificate
    from tools.stability.ou3_theorem.planar_complete_word_storage import certificate as word_storage_certificate
    from tools.stability.ou3_theorem.planar_innovation_storage import certificate as innovation_storage_certificate
    from tools.stability.ou3_theorem.information_shear_word import certificate as shear_word_certificate
    from tools.stability.ou3_theorem.measurement_frame import certificate as frame_loss_certificate
    for name, generate in (
        ("measurement-frame-loss-certificate.json",frame_loss_certificate),
        ("information-shear-word-certificate.json",shear_word_certificate),
        ("planar-innovation-storage-certificate.json",innovation_storage_certificate),
        ("planar-frontend-domain-certificate.json",frontend_domain_certificate),
        ("planar-complete-word-storage-certificate.json",word_storage_certificate),
        ("mahony-raw-normalization-certificate.json",raw_norm_certificate),
        ("planar-linked-riccati-mean-certificate.json",linked_calculus_certificate),
        ("planar-causal-calculus-certificate.json",causal_calculus_certificate),
        ("planar-mean-covariance-ports.json",mean_ports_certificate),
        ("planar-service-cell-certificate.json",planar_cell_certificate),
        ("planar-service-guard-certificate.json",planar_guard_certificate),
        ("planar-service-frontend-binding.json",planar_frontend_binding_certificate),
        ("planar-service-mahony-tube.json",planar_mahony_tube_certificate),
        ("planar-generated-port-structure.json",planar_generated_port_certificate),
        ("planar-joint-phase-cell-threshold.json",planar_joint_threshold_certificate),
        ("moving-quiet-compatibility-certificate.json",compatibility_certificate),
        ("imu-two-timescale-certificate.json",imu_certificate),
        ("word-diameter-certificate.json",diameter_certificate),
        ("world-frame-certificate.json",world_certificate),
        ("aw-covariance-ceiling-certificate.json",aw_ceiling_certificate),
        ("aw-tracking-certificate.json",aw_tracking_certificate),
        ("signed-injection-certificate.json",injection_certificate),
        ("aggregate-floor-certificate.json",aggregate_certificate),
        ("stationary-covariance-certificate.json",stationary_certificate),
        ("regime-continuation-feasibility.json",regime_diagnostic),
        ("regime-certificate.json",regime_certificate),
        ("moving-pivots-certificate.json",moving_certificate),
        ("gyro-bias-projection.json",gyro_certificate),
        ("ag-readout-certificate.json",readout_certificate),
        ("signed-temporal-certificate.json",signed_certificate),
        ("corrected-word-certificate.json",corrected_certificate),
        ("sampling-fidelity.json",sampling_certificate),
        ("lin-matrix-certificate.json",matrix_certificate),
        ("word-energy-audit.json",restricted_service_counterexample),
        ("block-factor-status.json",block_certificate),
        ("root-covariance-certificate.json",root_certificate),
        ("factor-nuisance-audit.json",audit_certificate),
        ("nuisance-upper-certificate.json",nuisance_certificate),
        ("sampled-capture-obstruction.json",witness_certificate),
    ):
        artifact=STATUS.parent/name
        if not artifact.is_file() or json.loads(artifact.read_text())!=generate():
            failures.append(f"committed {name} differs from exact reproduction")
    from tools.stability.ou3_theorem.planar_service_verify import verify as verify_planar
    try:
        verify_planar(*(json.loads((STATUS.parent/name).read_text()) for name in (
            "planar-service-stream-diagnostic.json", "planar-service-operation-audit.json",
            "planar-service-frechet-diagnostic.json")))
    except (OSError, ValueError, KeyError, TypeError) as error:
        failures.append("finite planar evidence verification failed: "+str(error))
    from tools.stability.ou3_theorem.planar_native_secants import verify_report as verify_native_secants
    from tools.stability.ou3_theorem.planar_moving_center import verify_report as verify_moving_center
    try:
        verify_moving_center(json.loads((STATUS.parent/'planar-moving-center-diagnostic.json').read_text()),
                            json.loads((STATUS.parent/'planar-service-stream-diagnostic.json').read_text())['stream_sha256'])
        verify_native_secants(json.loads((STATUS.parent/"planar-native-secants-diagnostic.json").read_text()))
        quotient=json.loads((STATUS.parent/"planar-quotient-mean-diagnostic.json").read_text())
        if quotient.get("result_type")!="FINITE DIAGNOSTIC ONLY" or quotient.get("quotient_dimension")!=20:
            raise ValueError("invalid finite physical quotient diagnostic")
        if quotient.get("stream_sha256")!=json.loads((STATUS.parent/"planar-service-stream-diagnostic.json").read_text())["stream_sha256"]:
            raise ValueError("quotient and native word are not linked")
        for key in ("complete_nonlinear_mean_derivative_computed","source_uniform_quotient_action_certified",
                    "joint_cell_forward_invariant","all_time_magnetic_service_verified","theorem_closed"):
            if quotient.get(key) is not False:raise ValueError("finite quotient promotion: "+key)
    except (OSError, ValueError, KeyError, TypeError) as error:
        failures.append("finite quotient/secant evidence verification failed: "+str(error))
    from tools.stability.ou3_theorem.construction_history_diagnostic import driver_source, zero_true_bias_storage_audit
    construction=json.loads((STATUS.parent/"construction-history-feasibility.json").read_text())
    if construction.get("generated_driver_sha256") != hashlib.sha256(driver_source().encode()).hexdigest():
        failures.append("construction native driver fingerprint changed")
    if construction.get("terminal_storage_certificate") != zero_true_bias_storage_audit(construction["native"]):
        failures.append("construction endpoint storage certificate differs from exact reproduction")
    from tools.stability.ou3_theorem.construction_mean_action import observer_source, instrument, HEADER, verify_summary
    mean=json.loads((STATUS.parent/"construction-mean-action.json").read_text())
    if mean.get("observer_sha256") != hashlib.sha256(observer_source().encode()).hexdigest():
        failures.append("construction mean observer fingerprint changed")
    if mean.get("instrumented_header_sha256") != hashlib.sha256(instrument((REPO/HEADER).read_text()).encode()).hexdigest():
        failures.append("construction mean header fingerprint changed")
    if mean.get("exact_summary") != verify_summary(mean["enclosure"], mean["committed_field"]):
        failures.append("construction mean action differs from rational factor reproduction")
    from tools.stability.ou3_theorem.signed_temporal_diagnostic import verify_diagnostic
    try:
        verify_diagnostic(json.loads((STATUS.parent/"signed-adjoint-diagnostic.json").read_text()))
    except (ValueError, KeyError, OSError) as error:
        failures.append("signed adjoint diagnostic verification failed: "+str(error))
    from tools.stability.ou3_theorem.signed_temporal_balance import verify_diagnostic as verify_balance
    try:
        verify_balance(json.loads((STATUS.parent/"signed-balance-diagnostic.json").read_text()))
    except (ValueError, KeyError, OSError) as error:
        failures.append("signed balance diagnostic verification failed: "+str(error))
    from tools.stability.ou3_theorem.moving_transport_source_diagnostic import verify_diagnostic as verify_groups
    try:
        verify_groups(json.loads((STATUS.parent/"moving-transport-source-feasibility.json").read_text()))
    except (ValueError, KeyError, OSError) as error:
        failures.append("same-cell diagnostic verification failed: "+str(error))
    from tools.stability.ou3_theorem.world_frame_source_diagnostic import verify_diagnostic as verify_world
    try:
        verify_world(json.loads((STATUS.parent/"world-frame-source-feasibility.json").read_text()))
    except (ValueError, KeyError, OSError) as error:
        failures.append("world-frame diagnostic verification failed: "+str(error))
    from tools.stability.ou3_theorem.aw_tracking_source_diagnostic import verify_diagnostic as verify_aw
    try:
        verify_aw(json.loads((STATUS.parent/"aw-tracking-source-feasibility.json").read_text()))
    except (ValueError, KeyError, OSError) as error:
        failures.append("AW tracking diagnostic verification failed: "+str(error))
    from tools.stability.ou3_theorem.information_ratio_source_diagnostic import verify_diagnostic as verify_ratio
    try:
        verify_ratio(json.loads((STATUS.parent/"information-ratio-source-feasibility.json").read_text()))
    except (ValueError, KeyError, OSError) as error:
        failures.append("information-ratio diagnostic verification failed: "+str(error))
    return {"validation_pass":not failures,"failures":failures,"base_main_commit":provenance["base_main_commit"],
            "shipping_behavior_authority":"source implementation","theorem_closed":expected["theorem_closed"]}

def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--output",type=Path); args=ap.parse_args()
    report=validate()
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(report,sort_keys=True)); return 0 if report["validation_pass"] else 1

if __name__=="__main__": raise SystemExit(main())
