"""The source-cover scalar checks must track the canonical physical envelope."""
import json
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/stability'))
import ou3_brmm_finite_window_primitive_qualification as PRIMITIVE
import ou3_p4_complete_brmm_source_cover_contract as COVER


class SourceCoverEnvelopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.envelope = json.loads(
            COVER.DEFAULT_DOMAIN.read_text())['complete_brmm_physical_envelope']
        cls.canonical = {
            'BRMM_uniform_V_m_mps': float(cls.envelope['wave_velocity_norm_upper_mps']),
            'BRMM_uniform_P_m_m': float(cls.envelope['wave_position_norm_upper_m']),
            'BRMM_DeltaS_3s_norm_upper_m_s': math.nextafter(
                PRIMITIVE.HORIZON_S * float(cls.envelope['wave_position_norm_upper_m']),
                math.inf),
        }

    def test_bounds_come_from_the_envelope_not_from_literals(self):
        self.assertEqual(COVER.V_M_UPPER_MPS, self.canonical['BRMM_uniform_V_m_mps'])
        self.assertEqual(COVER.P_M_UPPER_M, self.canonical['BRMM_uniform_P_m_m'])
        self.assertEqual(COVER.DELTA_S_3S_UPPER_M_S,
                         self.canonical['BRMM_DeltaS_3s_norm_upper_m_s'])

    def test_canonical_primitives_are_accepted(self):
        scalar = [f for f in COVER.validate(self.canonical) if f.endswith(' invalid')]
        self.assertEqual([], scalar)

    def test_nonpositive_or_oversized_primitives_are_rejected(self):
        for key, name in (('BRMM_uniform_V_m_mps', 'V_m'),
                          ('BRMM_uniform_P_m_m', 'P_m'),
                          ('BRMM_DeltaS_3s_norm_upper_m_s', 'DeltaS')):
            canonical = self.canonical[key]
            for value in (0.0, -canonical, math.nextafter(canonical, math.inf),
                          2.0 * canonical, float('nan'), float('inf')):
                with self.subTest(key=key, value=value):
                    forged = dict(self.canonical)
                    forged[key] = value
                    self.assertIn(name + ' invalid', COVER.validate(forged))


if __name__ == '__main__':
    unittest.main()
