"""Live magnetic sigma sensitivity on uninterrupted paired navigation replays."""
import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import t

from marine_mag_study import OUT, ROOT, parse
from sim_dataset import REFERENCE_CSV_SHA256_BY_CASE, digest

DRIVER = ROOT / "tests/kalman_tfg_shared_modes/magnetic_study/build/magnetic-heading-study"
KEY = ["case", "height", "seed", "scenario", "dip"]
METRICS = ["yaw_rms_deg", "abs_signed_yaw_deg", "centered_deg", "roll_rms_deg", "pitch_rms_deg", "disp_3d_rms_m"]


def simulate(task):
    case, height, seed, scenario, dip, factor, path = task
    env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "OCEAN_IMU_SIM_CASE", "MARINE_MAG_"))}
    env.update(OCEAN_IMU_SIM_CASE=case, W3D_SEED=str(seed), MARINE_MAG_SIGMA_FACTOR=str(factor))
    out = subprocess.check_output([DRIVER, path, "0", scenario, str(dip)], env=env, text=True)
    scores, _ = parse(out)
    return dict(zip(KEY, task[:5]), sigma_factor=factor, abs_signed_yaw_deg=abs(scores["yaw_mean_deg"]), **scores)


def panel(stage, factors, jobs):
    training = stage == "development"
    spectrum = "jonswap" if training else "pmstokes"
    cells = ([(0.27,"clean",66.5),(1.5,"hard",66.5),(4,"combined",66.5)] if training else
             [(0.27,"clean",45),(1.5,"combined",66.5),(4,"combined",55),(8.5,"combined",76)])
    seeds = (13,29) if training else (421,829,1609,3517)
    tasks = []
    for case in ("following","head"):
        for height, scenario, dip in cells:
            path, = (ROOT / ".sim-data/cases" / case).glob(f"wave_data_{spectrum}_H{height:.3f}_*.csv")
            assert digest(path) == REFERENCE_CSV_SHA256_BY_CASE[case][path.name]
            tasks.extend((case,height,seed,scenario,dip,factor,path) for seed in seeds for factor in factors)
    with ThreadPoolExecutor(jobs) as pool:
        frame = pd.DataFrame(pool.map(simulate,tasks))
    frame.to_csv(OUT / f"yaw-weight-{stage}.csv",index=False)
    return frame


def paired_scores(frame):
    original = frame[frame.sigma_factor==1].set_index(KEY)
    rows = []
    for factor, arm in frame.groupby("sigma_factor"):
        arm = arm.set_index(KEY).loc[original.index]
        for metric in METRICS:
            delta = arm[metric]-original[metric]
            draws = delta.groupby(level="seed").mean()
            margin = t.ppf(.975,len(draws)-1)*draws.std()/np.sqrt(len(draws))
            rows.append(dict(sigma_factor=factor,metric=metric,mean=arm[metric].mean(),
                delta=delta.mean(),pct=100*delta.mean()/original[metric].mean(),
                low=draws.mean()-margin,high=draws.mean()+margin,improved=int((delta < -1e-6).sum()),
                worse=int((delta > 1e-6).sum()),worst=delta.max(),
                worst_pct=(100*delta/original[metric]).max()))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs",type=int,default=4)
    args = parser.parse_args()
    development = panel("development",[.5,1.,2.,4.],args.jobs)
    base = pd.read_csv(OUT/"development.csv").query("variant==0").set_index(KEY)
    control = development.query("sigma_factor==1").set_index(KEY).loc[base.index]
    np.testing.assert_allclose(control.yaw_rms_deg,base.yaw_rms_deg,rtol=1e-6,atol=1e-6)
    train = paired_scores(development)
    by_factor = train.pivot(index="sigma_factor",columns="metric",values="delta")
    # Lock the ordinary-yaw candidate on development. A failed candidate is
    # still evaluated on new draws so a negative selection remains testable.
    eligible = []
    for factor in (.5,2.,4.):
        s = train[train.sigma_factor==factor].set_index("metric")
        if (by_factor.loc[factor,"yaw_rms_deg"]<0 and by_factor.loc[factor,"abs_signed_yaw_deg"]<0
                and s.at["disp_3d_rms_m","pct"]<=2 and s.at["disp_3d_rms_m","worst_pct"]<=10
                and max(s.at[k,"delta"] for k in ("roll_rms_deg","pitch_rms_deg"))<=.01
                and max(s.at[k,"worst"] for k in ("roll_rms_deg","pitch_rms_deg"))<=.1):
            eligible.append(factor)
    chosen = min(eligible or [.5,2.,4.], key=lambda k:by_factor.loc[k,"yaw_rms_deg"])
    print(f"Development qualified {eligible}, testing {chosen} on four fresh draws",flush=True)
    validation = panel("validation",[1.,chosen],args.jobs)
    stats = pd.concat([train.assign(stage="development"),paired_scores(validation).assign(stage="validation")])
    stats.to_csv(OUT/"yaw-weight-paired.csv",index=False)


if __name__ == "__main__":
    main()
