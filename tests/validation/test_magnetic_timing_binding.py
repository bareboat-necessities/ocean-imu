"""Keep legacy finite evidence distinct from the unqualified timed path."""
import hashlib
import json
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.magnetic_timing_legacy_binding import audit
from tools.stability.ou3_theorem.theorem_status import status_report

class MagneticTimingBindingTests(unittest.TestCase):
    def test_legacy_operation_identity_and_payload_preservation(self):
        result=audit()
        root=ROOT/'reports/results/ou3_stability'
        receipt=json.loads((root/'magnetic-timing-legacy-binding.json').read_text())
        self.assertEqual(result['sources'],receipt['sources'])
        self.assertFalse(receipt['timed_api_equivalence_claimed'])
        self.assertFalse(receipt['theorem_closed'])
        stream=json.loads((root/'planar-service-stream-diagnostic.json').read_text())
        self.assertEqual(stream['stream_sha256'],receipt['bitwise_native_control']['stream_sha256'])
        for name,bindings in receipt['rebound_records'].items():
            data=json.loads((root/name).read_text())
            self.assertEqual({k:data[k] for k in bindings['new_bindings']},bindings['new_bindings'])
            payload={k:v for k,v in data.items() if k not in bindings['new_bindings']}
            digest=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
            self.assertEqual(digest,bindings['unchanged_payload_sha256'],name)
    def test_timed_path_premises_remain_open(self):
        t=status_report()['magnetic_timing_preprocessing']
        self.assertFalse(t['physical_BMM150_conversion_timestamp_known'])
        self.assertFalse(t['legacy_reference_port_result_covers_endogenous_transport'])
        self.assertFalse(t['transport_state_noise_correlation_uniformly_bounded'])
        self.assertFalse(t['source_uniform_postprocessing_service_verified'])
        self.assertFalse(t['disturbance_residual_gating_enabled'])
        self.assertFalse(t['duplicate_host_observations_count_as_service'])
