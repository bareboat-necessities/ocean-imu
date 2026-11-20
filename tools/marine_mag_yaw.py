"""Magnetic weight and independent gyro-bias prior studies on uninterrupted replays."""
import argparse
import os
import subprocess
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import t

from marine_mag_study import OUT, ROOT, parse, physical_scenario, special
from sim_dataset import REFERENCE_CSV_SHA256_BY_CASE, digest

DRIVER = ROOT / "tests/kalman_tfg_shared_modes/magnetic_study/build/magnetic-heading-study"
KEY = ["case", "height", "seed", "scenario", "dip"]
METRICS = ["yaw_rms_deg", "abs_signed_yaw_deg", "centered_deg", "roll_rms_deg", "pitch_rms_deg", "disp_3d_rms_m"]


def simulate(task):
    case, height, seed, scenario, dip, factor, path, parameter = task
    env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "OCEAN_IMU_SIM_CASE", "MARINE_MAG_"))}
    env.update(OCEAN_IMU_SIM_CASE=case, W3D_SEED=str(seed))
    env["MARINE_MAG_GYRO_BIAS_FACTOR" if parameter == "gyro_bias_factor" else "MARINE_MAG_SIGMA_FACTOR"] = str(factor)
    fixture = physical_scenario(scenario, env)
    out = subprocess.check_output([DRIVER, path, "0", fixture, str(dip)], env=env, text=True)
    scores, _ = parse(out)
    row = dict(zip(KEY, task[:5]))
    row[parameter] = factor
    return dict(row, abs_signed_yaw_deg=abs(scores["yaw_mean_deg"]), **scores)


def panel(stage, factors, jobs, parameter, prefix, validation_seeds):
    training = stage == "development"
    spectrum = "jonswap" if training else "pmstokes"
    cells = ([(0.27,"clean",66.5),(1.5,"hard",66.5),(4,"combined",66.5)] if training else
             [(0.27,"clean",45),(1.5,"combined",66.5),(4,"combined",55),(8.5,"combined",76)])
    seeds = (13,29) if training else validation_seeds
    tasks = []
    for case in ("following","head"):
        for height, scenario, dip in cells:
            path, = (ROOT / ".sim-data/cases" / case).glob(f"wave_data_{spectrum}_H{height:.3f}_*.csv")
            assert digest(path) == REFERENCE_CSV_SHA256_BY_CASE[case][path.name]
            tasks.extend((case,height,seed,scenario,dip,factor,path,parameter) for seed in seeds for factor in factors)
    with ThreadPoolExecutor(jobs) as pool:
        frame = pd.DataFrame(pool.map(simulate,tasks))
    frame.to_csv(OUT / f"{prefix}-{stage}.csv",index=False)
    return frame


def gyro_stress(jobs):
    fixtures = pd.read_csv(OUT/"stress-inputs.csv")
    tasks = []
    with tempfile.TemporaryDirectory() as temporary:
        for row in fixtures.itertuples(index=False):
            data = ROOT/".sim-data/cases"/row.case
            scenario = row.scenario
            path = special(data,row.case,Path(temporary),scenario) if scenario in ("turn","weak","roundtrip") else data/row.file
            if digest(path) != row.sha256:
                raise ValueError(f"changed stress input {path}")
            if scenario in ("turn","weak") and digest(path.with_suffix(".track")) != row.track_sha256:
                raise ValueError(f"changed scoring track {path}")
            tasks.extend((row.case,1.5,3037,scenario,66.5,k,path,"gyro_bias_factor") for k in (1.,4.))
        with ThreadPoolExecutor(jobs) as pool:
            scored = pd.DataFrame(pool.map(simulate,tasks))
    native = scored[scored.gyro_bias_factor==1].set_index(KEY)
    reference = pd.read_csv(OUT/"stress.csv").query("variant==0").set_index(KEY).loc[native.index]
    np.testing.assert_allclose(native.yaw_rms_deg,reference.yaw_rms_deg,rtol=1e-6,atol=1e-6)
    scored.to_csv(OUT/"gyro-bias-stress.csv",index=False)



def paired_scores(frame, parameter):
    original = frame[frame[parameter]==1].set_index(KEY)
    rows = []
    for factor, arm in frame.groupby(parameter):
        arm = arm.set_index(KEY).loc[original.index]
        for metric in METRICS:
            delta = arm[metric]-original[metric]
            draws = delta.groupby(level="seed").mean()
            margin = t.ppf(.975,len(draws)-1)*draws.std()/np.sqrt(len(draws))
            rows.append(dict({parameter:factor},metric=metric,mean=arm[metric].mean(),
                delta=delta.mean(),pct=100*delta.mean()/original[metric].mean(),
                low=draws.mean()-margin,high=draws.mean()+margin,improved=int((delta < -1e-6).sum()),
                worse=int((delta > 1e-6).sum()),worst=delta.max(),
                worst_pct=(100*delta/original[metric]).max()))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs",type=int,default=4)
    parser.add_argument("--gyro-bias",action="store_true",help="scale learner gyro-bias prior and random walk")
    parser.add_argument("--replicate-gyro-bias",action="store_true",help="independently replicate the locked fourfold gyro-bias factor")
    parser.add_argument("--gyro-stress",action="store_true",help="replay the fixed stress conditions with the locked factor")
    args = parser.parse_args()
    if args.gyro_stress:
        gyro_stress(args.jobs)
        return
    if args.replicate_gyro_bias:
        seeds = (5003,5011,5021,5023,5039,5051,5059,5077,5081,5087,5099,5101)
        validation = panel("replication",[1.,4.],args.jobs,"gyro_bias_factor","gyro-bias",seeds)
        stats = paired_scores(validation,"gyro_bias_factor")
        stats.to_csv(OUT/"gyro-bias-replication-paired.csv",index=False)
        return
    if args.gyro_bias:
        parameter, prefix = "gyro_bias_factor", "gyro-bias"
        alternatives, validation_seeds = (0.,.25,4.), (547,1151,2371,4099)
    else:
        parameter, prefix = "sigma_factor", "yaw-weight"
        alternatives, validation_seeds = (.5,2.,4.), (421,829,1609,3517)
    factors = sorted([1., *alternatives])
    development = panel("development",factors,args.jobs,parameter,prefix,validation_seeds)
    base = pd.read_csv(OUT/"development.csv").query("variant==0").set_index(KEY)
    control = development[development[parameter]==1].set_index(KEY).loc[base.index]
    np.testing.assert_allclose(control.yaw_rms_deg,base.yaw_rms_deg,rtol=1e-6,atol=1e-6)
    train = paired_scores(development,parameter)
    by_factor = train.pivot(index=parameter,columns="metric",values="delta")
    # Lock the ordinary-yaw candidate on development. A failed candidate is
    # still evaluated on new draws so a negative selection remains testable.
    eligible = []
    for factor in alternatives:
        s = train[train[parameter]==factor].set_index("metric")
        if (by_factor.loc[factor,"yaw_rms_deg"]<0 and by_factor.loc[factor,"abs_signed_yaw_deg"]<0
                and s.at["disp_3d_rms_m","pct"]<=2 and s.at["disp_3d_rms_m","worst_pct"]<=10
                and max(s.at[k,"delta"] for k in ("roll_rms_deg","pitch_rms_deg"))<=.01
                and max(s.at[k,"worst"] for k in ("roll_rms_deg","pitch_rms_deg"))<=.1):
            eligible.append(factor)
    chosen = min(eligible or alternatives, key=lambda k:by_factor.loc[k,"yaw_rms_deg"])
    print(f"Development qualified {eligible}, testing {chosen} on four fresh draws",flush=True)
    validation = panel("validation",[1.,chosen],args.jobs,parameter,prefix,validation_seeds)
    stats = pd.concat([train.assign(stage="development"),paired_scores(validation,parameter).assign(stage="validation")])
    stats.to_csv(OUT/f"{prefix}-paired.csv",index=False)


if __name__ == "__main__":
    main()
