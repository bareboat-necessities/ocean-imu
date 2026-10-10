"""Mechanical legacy-API operation audit for magnetic preprocessing.

Partial-evaluate ONLY the original updateMag entry (null observation/covariance
and the unchanged magnetic-correction callback). Remove read-only diagnostics
and new, uncalled overloads/accessors. The complete resulting three headers
must match the pre-change source token-for-token. This does not cover the new
timed entry, compiler roundoff in a different profile, or hardware sampling.
"""
from pathlib import Path
import hashlib
import re

ROOT = Path(__file__).resolve().parents[3]
BASE_COMMIT = '2c534feb423a627c40f475daed40b4274d3dd2b8'
PATHS = ('src/kalman_ou_iii/Kalman3D_Wave_OU_III.h',
         'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h',
         'src/kalman_common/ProxyStartupFusion.h')
EXPECTED = {'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h': 'f87f7b5d533701df304744b36b5e08005dc4688e78e60b9e1ffb11b3a78cbe2a', 'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h': 'c797f8f614f38e616ce095dc4bbbcab11c4aa724b176f0a67f5e524f0fa2af76', 'src/kalman_common/ProxyStartupFusion.h': 'f1b7878da196a2439d81e0343414456c5398ad27f6365d0235c8ea45415377a1'}


def tokens(source):
    source = re.sub(r'/\*.*?\*/|//[^\n]*', '', source, flags=re.S)
    return re.sub(r'\s+', '', source)


def project(path, source):
    if path == PATHS[0]:
        start = source.index('    // Full covariance in physical body axes')
        end = source.index('    // Extended-only API:', start)
        source = source[:start] + '\n' + source[end:]
        start = source.index('    [[nodiscard]] Vector3 gyroscope_bias_body()')
        end = source.index('    [[nodiscard]] Vector3 get_acc_bias()', start)
        source = source[:start] + source[end:]
    elif path == PATHS[1]:
        if 'void updateMag(const Eigen::Vector3f& mag_body_ned) {\n        updateMagImpl_(mag_body_ned, nullptr);\n    }' not in source:
            raise ValueError('legacy magnetic covariance selector changed')
        start = source.index('    void updateMag(const Eigen::Vector3f& mag_body_ned) {')
        end = source.index('        if (!with_mag_ || !mekf_) return;', start)
        source = source[:start] + '    void updateMag(const Eigen::Vector3f& mag_body_ned) {\n' + source[end:]
        source = source.replace('        if (covariance_body) mekf_->measurement_update_mag_only(mag_body_ned, *covariance_body);\n        else ', '        ')
        source = source.replace('public:\n    void setWithMag', '    void setWithMag')
    elif path == PATHS[2]:
        # Verify the actual legacy entry really selects null observations and
        # exactly the original correction. The timed callback is out of scope.
        required = '''updateMagImpl_(mag_body_ned, nullptr, [this](const Eigen::Vector3f& field) {
            impl_.updateMag(field);
            return impl_.mekf().lastMagDiag().accepted;
        });'''
        if required not in source:
            raise ValueError('legacy magnetic callback changed')
        source = source.replace('#include "util/MagneticRotation.h"', '')
        start = source.index('    void updateMag(const Eigen::Vector3f& mag_body_ned) {')
        end = source.index('        if (!begun_ || !cfg_.with_mag) return;', start)
        source = source[:start] + '    void updateMag(const Eigen::Vector3f& mag_body_ned) {\n' + source[end:]
        source = source.replace('        const float sample_t = observation ? t_ - observation->rotation.seconds : t_;\n', '')
        source = source.replace(', observation, sample_t)', ')')
        source = source.replace('observation ? observation->attitude_bw : attitudeReferenceQuat_()', 'attitudeReferenceQuat_()')
        source = source.replace('observation ? observation->accel : last_acc_body_ned_', 'last_acc_body_ned_')
        source = source.replace('observation ? observation->gyro : last_gyro_body_ned_', 'last_gyro_body_ned_')
        source = re.sub(r'observation \? yawRemovedBoatQuat\(observation->proxy_bw\)\s*:\s*impl_\.startupProxyTiltQuat\(\)',
                        'impl_.startupProxyTiltQuat()', source)
        source = source.replace('last_mag_applied_ = correct(', 'impl_.updateMag(')
        source = source.replace('public:\n    bool hasMagNorthLock', '    bool hasMagNorthLock')
        source = source.replace('        last_mag_applied_ = false;\n', '')
        source = source.replace('    bool last_mag_applied_ = false;\n', '')
        source = re.sub(r'const Eigen::Vector3f& mag_body_ned,\s*const ocean_imu::magnetic::Observation\* observation, float sample_t',
                        'const Eigen::Vector3f& mag_body_ned', source)
        source = re.sub(r'\bsample_t\b', 't_', source)
    return source


def audit():
    result = {}
    for path in PATHS:
        source = (ROOT/path).read_text()
        sha = hashlib.sha256(tokens(project(path, source)).encode()).hexdigest()
        if sha != EXPECTED[path]:
            raise ValueError('legacy operation projection changed: ' + path)
        result[path] = {'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
                        'legacy_operation_tokens_sha256': sha}
    return {'baseline_commit': BASE_COMMIT, 'legacy_api_operation_projection_verified': True,
            'timed_api_equivalence_claimed': False, 'theorem_closed': False, 'sources': result}
