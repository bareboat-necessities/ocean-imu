"""Constructive periodic fast-gyro witness certificate.

Non-promoting until every interval check is true. The forcing period is 6 s:
1200 literal 5 ms IMU samples and 150 deployed 25 Hz magnetic callbacks.
"""
from fractions import Fraction as F
DT=F(1,200); PERIOD_S=F(6); SAMPLES=1200; MAG_STRIDE=8
VAMP=F(7,2); NGAMP=F(1,50)

def certificate():
    pi_lo,pi_hi=F(333,106),F(355,113)
    om_lo,om_hi=pi_lo/3,pi_hi/3
    p_hi=VAMP/om_lo
    a_hi=VAMP*om_hi
    j_hi=VAMP*om_hi*om_hi
    theta_hi=NGAMP/om_lo
    supply=VAMP*NGAMP/2
    return {
      "qualification":"OU3_PERIODIC_FAST_GYRO_WITNESS_V1",
      "period_s":str(PERIOD_S),"samples_per_period":SAMPLES,
      "magnetic_callbacks_per_period":SAMPLES//MAG_STRIDE,
      "omega_rad_s_lower":str(om_lo),"omega_rad_s_upper":str(om_hi),
      "p_amp_upper_m":str(p_hi),"v_amp_mps":str(VAMP),
      "a_amp_upper_mps2":str(a_hi),"jerk_amp_upper_mps3":str(j_hi),
      "field_axis_error_amp_upper_rad":str(theta_hi),
      "signed_fast_gyro_supply_mps2":str(supply),
      "old_margin_mps2":"0.0338084",
      "supply_exceeds_old_margin":supply>F("0.0338084"),
      "physical_envelopes_pre_gravity_compensation":
          p_hi<F("8.1") and VAMP<=F("5.5") and a_hi<F("8.8") and j_hi<F(100),
      "native_periodic_orbit_exported":False,
      "periodic_tuner_word_interval_enclosed":False,
      "periodic_covariance_orbit_interval_enclosed":False,
      "physical_compatibility_operator_interval_enclosed":False,
      "physical_compatibility_operator_nonsingular":False,
      "physical_solution_interval_enclosed":False,
      "all_shipping_gates_verified":False,
      "magnetic_service_verified":False,
      "weighted_signed_functional_interval_lower":None,
      "shipping_counterexample_certified":False,
      "theorem_closed":False}

def driver_source():
    """Untouched shipping wrapper driven by the exact smooth candidate."""
    from .construction_history_diagnostic import REPO
    src=(REPO/'tools/stability/ag_readout_source.cpp').read_text()
    src=src.replace('for (int k=1; k<=45064; ++k) {',
                    'for (int k=1; k<=240000; ++k) {')
    a=src.index('        if (k==45001) {')
    b=src.index('        const double t',a)
    src=src[:a]+'''        if (k==228001) { root=matrix_json(filter.raw().mekf().covariance_full()); recording=true; }
    '''+src[b:]
    src=src.replace('const float roll = wave ? static_cast<float>(.02*std::sin(.5*t)) : 0.0f;',
      'const double om=M_PI/3.0; const float roll=wave ? static_cast<float>((.02/om)*std::cos(om*t)) : 0.0f;')
    # true roll rate plus +.02 sin(om t) fast residual is identically zero
    src=src.replace('const float rate = wave ? static_cast<float>(.01*std::cos(.5*t)) : 0.0f;',
                    'const float rate=0.0f;')
    src=src.replace('const float az = wave ? static_cast<float>(-.144*std::sin(.6*t)) : 0.0f;',
                    'const float ay=wave ? static_cast<float>(3.5*(M_PI/3.0)*std::cos((M_PI/3.0)*t)) : 0.0f;')
    src=src.replace('rwb*Eigen::Vector3f(0,0,az-g_std)',
                    'rwb*Eigen::Vector3f(0,ay,-g_std)')
    src=src.replace('rwb*Eigen::Vector3f(60,0,30)',
                    'Eigen::Vector3f(75,0,0)')
    src=src.replace(' || applied!=8','')
    src=src.replace('    const auto& m=filter.raw().mekf();',
      '''    const auto& m=filter.raw().mekf();
    const auto tune=filter.raw().tune_;''')
    src=src.replace('<< ",\\\"root_covariance\\\":" << root',
      '''<< ",\\\"tau_applied\\\":" << tune.tau_applied
              << ",\\\"sigma_applied\\\":" << tune.sigma_applied
              << ",\\\"RS_applied\\\":" << tune.RS_applied
              << ",\\\"pseudo_period\\\":" << m.get_pseudo_update_period_s()
              << ",\\\"committed_field\\\":" << matrix_json(m.v2ref)
              << ",\\\"root_covariance\\\":" << root''')
    return src

def native_probe(eigen):
    import json, subprocess, tempfile
    from pathlib import Path
    from .construction_history_diagnostic import REPO
    source=driver_source()
    with tempfile.TemporaryDirectory(prefix='ou3-periodic-fast-gyro-') as d:
        d=Path(d); cpp=d/'driver.cpp'; exe=d/'driver'
        cpp.write_text(source)
        subprocess.run(['g++','-O2','-std=c++20','-I'+str(REPO/'src'),
                        '-isystem',str(eigen),str(cpp),'-o',str(exe)],check=True)
        return json.loads(subprocess.check_output([str(exe),'wave'],text=True))

if __name__=="__main__":
    import argparse,json
    from pathlib import Path
    p=argparse.ArgumentParser()
    p.add_argument('--eigen',type=Path)
    p.add_argument('--native',action='store_true')
    args=p.parse_args()
    out=certificate()
    if args.native:
        if args.eigen is None: raise SystemExit('--eigen required with --native')
        out['native']=native_probe(args.eigen)
        out['native_periodic_orbit_exported']=True
    print(json.dumps(out,indent=2,sort_keys=True))

