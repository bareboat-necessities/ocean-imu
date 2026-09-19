"""Fail-closed assembly of LIN factors and complete A21 word obligations.

A LIN endpoint factor and restricted magnetic information do not constitute a
full-state contraction certificate. No hard-coded AG/BA probe is promoted.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from tools.stability.ou3_theorem.lin_path_certificate import certificate as lin_certificate
from tools.stability.ou3_theorem.lin_matrix_certificate import certificate as matrix_certificate
from tools.stability.ou3_theorem.word_energy import restricted_service_counterexample


def certificate() -> dict:
    lin=lin_certificate()
    matrix=matrix_certificate()
    audit=restricted_service_counterexample()
    return {
        'qualification':'OU3_A21_BLOCK_FACTOR_ASSEMBLY_V2',
        'verified':False, 'mu_cov':None, 'rho0':None,
        'lin_path_verified':lin['verified'],
        'lin_matrix_verified':matrix['verified'],
        'ell_LIN':lin['factor_floor'],
        'raw_translation_normalization_lower':matrix['normalization_raw_translation_lower'],
        'full_word_energy_algebra_verified':audit['verified'],
        'restricted_service_lifting_valid':False,
        'constructive_full_A21_mu_rho_enclosure':False,
        'open_obligations':[
            'source-uniform full transported loss including cross-coordinate cancellation',
            'recurring AG and BA factors through literal corrections and resets',
            'complete finite-error projection/reset/tuner composition',
            'float32 supply, coercivity and every-prefix retention',
        ],
        'float32_whole_word_supply_closed':False, 'theorem_closed':False,
    }


def main() -> int:
    report=certificate()
    print(json.dumps(report,indent=2,sort_keys=True))
    # This command reproduces partial evidence. It is not a theorem gate.
    return 0 if report['lin_matrix_verified'] and report['full_word_energy_algebra_verified'] else 1


if __name__=='__main__':
    raise SystemExit(main())
