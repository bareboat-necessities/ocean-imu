"""Mechanical legacy-API operation audit for magnetic preprocessing.

Partial-evaluate the original updateMag entry (null observation/covariance and
the unchanged magnetic-correction callback) with attitude aiding and fresh
acceleration always enabled and no external reacquisition. Remove read-only
diagnostics and uncalled overloads/accessors. The complete resulting three
headers must match the pre-change source token-for-token. This does not cover
the motion guard, recovery, timed entry, different-profile compiler roundoff,
or hardware sampling.
"""
from pathlib import Path
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[3]
BASE_COMMIT = '2c534feb423a627c40f475daed40b4274d3dd2b8'
PATHS = ('src/kalman_ou_iii/Kalman3D_Wave_OU_III.h',
         'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h',
         'src/kalman_common/ProxyStartupFusion.h')
EXPECTED = {'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h': '95aa5cd90dd8d0e1345964095cc6d851d5a337d9b1ce2906fe29f5f1a93999e9', 'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h': '41a2b0524f0174880f9045f8e8128ebb45260e39b659be5d8c48354c59a69f0d', 'src/kalman_common/ProxyStartupFusion.h': 'b2b061bbfc22a4ebe4236412ae179909a991dffcd0460bc370749e96e1c8d1cd'}


_TOKEN = re.compile(
    r'(?P<ignored>\s+|//[^\n]*|/\*.*?\*/)'
    r'|(?:u8|u|U|L)?R"(?P<delimiter>[^ ()\\\t\r\n]{0,16})\(.*?\)(?P=delimiter)"'
    r'|(?:u8|u|U|L)?(?:"(?:\\.|[^"\\\r\n])*"|\'(?:\\.|[^\'\\\r\n])*\')'
    r'|(?:[0-9]|\.[0-9])(?:[eEpP][+-]|[a-zA-Z0-9_.\'])*'
    r'|[a-zA-Z_][a-zA-Z0-9_]*'
    r'|%:%:|>>=|<<=|<=>|->\*|\.\.\.|##|::|\.\*|->|\+\+|--|<<|>>'
    r'|<=|>=|==|!=|&&|\|\||[+*/%&|^!-]=|<:|:>|<%|%>|%:|.',
    re.S,
)


def tokens(source):
    """Serialize lexical tokens, retaining literal bytes and token boundaries.

    Comments/spacing outside literals are insignificant to this source audit.
    This is not preprocessing or a general C++ semantic-equivalence proof.
    """
    return json.dumps([match.group() for match in _TOKEN.finditer(source)
                       if match.group('ignored') is None], separators=(',', ':'))


def replace_once(source, old, new=''):
    if source.count(old) != 1:
        raise ValueError('ordinary-aiding projection anchor changed: ' + old[:80])
    return source.replace(old, new, 1)


def ordinary_aiding(path, source):
    """Specialize only the declared true flags; the final token hash is fixed."""
    if path == PATHS[0]:
        # The setters are excluded from this execution domain. Verify the
        # actual initializers before eliminating the constant true guards.
        for flag in ('attitude_aiding_enabled_', 'accelerometer_fresh_'):
            source = replace_once(source, '    bool ' + flag + ' = true;\n')
        start = source.index('    // Device motion policy:')
        end = source.index('    void time_update(', start)
        source = source[:start] + source[end:]
        source = replace_once(source, '''        if (!attitude_aiding_enabled_) {
            K.template topRows<3>().setZero();
            if constexpr (with_gyro_bias) K.template middleRows<3>(3).setZero();
            if constexpr (with_accel_bias) K.template middleRows<3>(OFF_BA).setZero();
        }
''')
        source = replace_once(source, '    if (!accelerometer_fresh_) return;\n')
        source = replace_once(source, '    if (!attitude_aiding_enabled_) return;\n')
        source = replace_once(source, '''        const bool held = !attitude_aiding_enabled_ &&
            (i < BASE_N || (with_accel_bias && i >= OFF_BA));
        K(i) = held ? T(0) : PCt(i) / S;''', '        K(i) = PCt(i) / S;')
    elif path == PATHS[1]:
        source = replace_once(source, '''    void setAttitudeAidingEnabled(bool enabled) {
        mekf_->set_attitude_aiding_enabled(enabled);
        this->vertical_accel_comp_.setAidingEnabled(enabled);
    }
    void setAccelerometerFresh(bool fresh) { mekf_->set_accelerometer_fresh(fresh); }
''')
        source = replace_once(source, 'drive_mekf && mekf_->attitude_aiding_enabled() &&', 'drive_mekf &&')
    elif path == PATHS[2]:
        source = replace_once(source, 'if (impl_.mekf().attitude_aiding_enabled()) gravity_gate_.step(', 'gravity_gate_.step(')
        source = replace_once(source, 'if (!begun_ || !impl_.mekf().attitude_aiding_enabled()) return;', 'if (!begun_) return;')
    return source


def project(path, source):
    source = ordinary_aiding(path, source)
    if path == PATHS[0]:
        source = source.replace('#include "util/MeasurementCovariance.h"', '')
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
            'execution_scope': 'legacy API; aiding always enabled; every acceleration fresh; no external reacquisition',
            'guarded_or_recovery_equivalence_claimed': False,
            'timed_api_equivalence_claimed': False, 'theorem_closed': False, 'sources': result}
