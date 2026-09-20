"""Carried construction stress history; finite evidence, never a uniform proof.

One fixed continuous physical history feeds the literal wrapper from begin().
No estimator setter, reseed, covariance replacement or artificial root is used.
The original observer driver is copied into a temporary directory; this probe
uses the untouched shipping headers and no observer instrumentation.
"""
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]

def driver_source():
    s=(REPO/'tools/stability/ag_readout_source.cpp').read_text()
    s=s.replace('for (int k=1; k<=45064; ++k) {', 'double tilt_sum=0, tilt_max=0, bg_max=0, aw_max=0, cross_min=1; int scored=0;\n    for (int k=1; k<=120000; ++k) {')
    a=s.index('        if (k==45001) {');b=s.index('        const double t',a)
    s=s[:a]+'''        if (k==45001) root=matrix_json(filter.raw().mekf().covariance_full());
    ''' +s[b:]
    s=s.replace('const float roll = wave ? static_cast<float>(.02*std::sin(.5*t)) : 0.0f;', 'const float roll=0.0f;')
    s=s.replace('const float rate = wave ? static_cast<float>(.01*std::cos(.5*t)) : 0.0f;', 'const float rate=0.0f;')
    s=s.replace('const float az = wave ? static_cast<float>(-.144*std::sin(.6*t)) : 0.0f;', 'const float az = wave ? static_cast<float>(6*std::sin(2*t)) : 0.0f;')
    s=s.replace('rwb*Eigen::Vector3f(0,0,az-g_std)', 'rwb*Eigen::Vector3f(az,0,az-g_std)')
    s=s.replace('rwb*Eigen::Vector3f(60,0,30)', 'rwb*Eigen::Vector3f(45,0,45)')
    s=s.replace('if (recording && filter.raw().mekf().lastMagDiag().accepted) ++applied;', 'if (filter.raw().mekf().lastMagDiag().accepted) ++applied;')
    s=s.replace('        if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;', '''        if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;
            const auto& m=filter.raw().mekf();
            bg_max=std::max(bg_max,static_cast<double>(m.xext.segment<3>(3).norm()));
            aw_max=std::max(aw_max,static_cast<double>(m.xext.segment<3>(15).norm()));
            if(k>80000) {
                const auto r=m.quaternion_boat().toRotationMatrix();
                const double tilt=std::acos(std::clamp(static_cast<double>(r(2,2)),-1.0,1.0))*180/M_PI;
                tilt_sum+=tilt; tilt_max=std::max(tilt_max,tilt); ++scored;
                const Eigen::Vector3f f=m.xext.segment<3>(15)-Eigen::Vector3f(0,0,g_std);
                const Eigen::Vector3f b=m.v2ref;
                cross_min=std::min(cross_min,static_cast<double>(f.cross(b).norm()/(f.norm()*b.norm())));
            }''')
    s=s.replace(' || applied!=8','')
    s=s.replace('    std::cout << "{\\\"live_step\\\":" << live', '''    std::cout << std::setprecision(17) << "{\\\"tail_tilt_mean_deg\\\":" << tilt_sum/scored
                  << ",\\\"tail_tilt_max_deg\\\":" << tilt_max << ",\\\"construction_bg_max\\\":" << bg_max
                  << ",\\\"construction_aw_max\\\":" << aw_max << ",\\\"tail_force_field_sin_min\\\":" << cross_min
                  << ",\\\"live_step\\\":" << live''')
    s=s.replace('static_cast<double>(k)*.005', 'static_cast<double>(k)*static_cast<double>(.005f)')
    for anchor in ('k<=120000', 'const float roll=0.0f;', 'const float rate=0.0f;',
                   '6*std::sin(2*t)', 'rwb*Eigen::Vector3f(az,0,az-g_std)',
                   'rwb*Eigen::Vector3f(45,0,45)', 'k>80000'):
        if s.count(anchor) != 1:
            raise ValueError('construction profile binding changed: '+anchor)
    if 'recording=true' in s:
        raise ValueError('construction probe must not alter or instrument the estimator')
    return s


def run(eigen):
    import hashlib
    import json
    import subprocess
    import tempfile
    from fractions import Fraction as F

    source = driver_source()
    with tempfile.TemporaryDirectory(prefix='ou3-construction-') as directory:
        driver, binary = Path(directory)/'construction.cpp', Path(directory)/'construction'
        driver.write_text(source)
        subprocess.run(['g++', '-O2', '-std=c++20', '-I'+str(REPO/'src'),
                        '-isystem', str(eigen), str(driver), '-o', str(binary)], check=True)
        observed = json.loads(subprocess.check_output([str(binary), 'wave'], text=True))
    # For p=-(3/2) sin(2t)(1,0,1): v=-3cos(2t)u, a=6sin(2t)u,
    # jerk=12cos(2t)u; a primitive is (3/4)cos(2t)u. These are all-time
    # continuous truth identities, independently of the finite native replay.
    squared = {'p': F(9, 2), 'v': F(18), 'a': F(72), 'jerk': F(288),
               'primitive_span': F(9, 2), 'omega': F(0)}
    limits = {'p': F('8.1'), 'v': F('5.5'), 'a': F('8.8'), 'jerk': F(100),
              'primitive_span': F(1100), 'omega': F('0.6108652381980153')}
    if not all(squared[k] <= limits[k]**2 for k in squared):
        raise ArithmeticError('physical history exceeds unchanged limits')
    return {'qualification': 'OU3_CARRIED_CONSTRUCTION_STRESS_DIAGNOSTIC_V1',
            'physical_history': 'p=-(3/2)sin(2t)(1,0,1), R=I, B=(45,0,45), zero true biases',
            'all_time_physical_squared_envelopes': {k: str(v) for k, v in squared.items()},
            'physical_envelopes_and_kinematics_verified_analytically': True,
            'physical_acceleration_DC': 'zero',
            'generated_driver_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'sampling': 'dt=.005f; physical time k*double(.005f); actual wrapper float clock retained',
            'tail_samples': 'steps 80001 through 120000',
            'native': observed,
            'all_time_magnetic_service_verified': False,
            'real_arithmetic_source_trajectory_enclosed': False,
            'six_degree_eventual_capture_refuted': False,
            'source_uniform_action_verified': False, 'theorem_closed': False}


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen', type=Path, default=Path('/usr/include/eigen3'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(args.eigen), indent=2, sort_keys=True)+'\n')
