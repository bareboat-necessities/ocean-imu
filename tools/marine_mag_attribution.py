"""Match physical inputs across magnetic installation, packet noise and learner settings."""
import argparse
import csv
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from statistics import mean, stdev


# Fixed before running. One sensor/installation/initialization draw is the
# independent unit. Each draw sees the same eight physical records in every arm.
DRAWS = (6011, 6029, 6043, 6067, 6089, 6101, 6121, 6133)
HEIGHT_DIP = ((.27, 45), (1.5, 66.5), (4., 55), (8.5, 76))
WEAK_DRAWS = (7001, 7013, 7027, 7043, 7057, 7079, 7103, 7121)
WEAK_HEIGHT_DIP = ((.27, 45), (.27, 76))
SCENARIOS = ("ideal", "clean", "combined")
VARIANTS = (-1, 0, 11)
METRICS = ("yaw_rms_deg", "abs_signed_yaw_deg", "centered_deg", "roll_rms_deg",
           "pitch_rms_deg", "disp_3d_rms_m")
CONTRASTS = {
    "screen_on_distorted": (("combined", 0), ("combined", 11)),
    "screen_on_clean": (("clean", 0), ("clean", 11)),
    "screen_on_ideal": (("ideal", 0), ("ideal", 11)),
    "learner_on_distorted": (("combined", -1), ("combined", 0)),
    "learner_on_clean": (("clean", -1), ("clean", 0)),
    "learner_on_ideal": (("ideal", -1), ("ideal", 0)),
    "installation_with_learner": (("clean", 0), ("combined", 0)),
    "installation_without_learner": (("clean", -1), ("combined", -1)),
    "packet_noise_without_learner": (("ideal", -1), ("clean", -1)),
    "packet_noise_with_learner": (("ideal", 0), ("clean", 0)),
}


def simulate(task):
    from marine_mag_study import BINARY, parse, physical_scenario

    case, height, dip, seed, scenario, variant, path = task
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("W3D_", "MARINE_MAG_", "OCEAN_IMU_SIM_CASE"))}
    env.update(W3D_SEED=str(seed), OCEAN_IMU_SIM_CASE=case)
    fixture = physical_scenario(scenario, env)
    output = subprocess.check_output([BINARY, path, str(variant), fixture, str(dip)],
                                     env=env, text=True)
    scores, _ = parse(output)
    events = [int(line.split()[1]) for line in output.splitlines() if line.startswith("INSTALLATION_RELEASE ")]
    if len(events) > 1 or (events and variant != 11):
        raise ValueError(f"Unexpected installation release: {task}, {events}")
    released_at = events[0] if events else -1

    scores["abs_signed_yaw_deg"] = abs(scores["yaw_mean_deg"])
    return dict(case=case, height=height, dip=dip, seed=seed, scenario=scenario,
                variant=variant, release_s=released_at if variant == 11 else None,
                **{name: scores[name] for name in (*METRICS, "yaw_mean_deg")})


def navigation_scores(row):
    scores = {name: float(row[name]) for name in (*METRICS, "yaw_mean_deg")}
    if (not all(map(math.isfinite, scores.values())) or any(scores[k] < 0 for k in METRICS)
            or abs(scores["yaw_mean_deg"]) != scores["abs_signed_yaw_deg"]
            or not math.isclose(scores["yaw_rms_deg"]**2, scores["centered_deg"]**2 + scores["yaw_mean_deg"]**2,
                                abs_tol=1e-5)):
        raise ValueError(f"Inconsistent navigation scores: {row}")
    if int(row["variant"]) == 11:
        release = float(row["release_s"])
        if release != -1 and not (15 <= release <= 1195 and release % 5 == 0):
            raise ValueError(f"Invalid release epoch: {row}")
    elif row["release_s"] not in (None, ""):
        raise ValueError(f"Unexpected release epoch for control: {row}")
    return scores


def summarize(rows, draw_ids, conditions):
    cells = {}
    expected = {(case, h, dip, seed) for case, (h, dip), seed
                in product(("following", "head"), conditions, draw_ids)}
    for row in rows:
        key = (row["case"], float(row["height"]), float(row["dip"]), int(row["seed"]))
        arm = row["scenario"], int(row["variant"])
        if key not in expected or arm in cells.setdefault(key, {}):
            raise ValueError(f"Unexpected or duplicate score: {key}, {arm}")
        cells[key][arm] = navigation_scores(row)
    if set(cells) != expected or any(set(arms) != set(product(SCENARIOS, VARIANTS)) for arms in cells.values()):
        raise ValueError("Incomplete installation screen")

    summaries = []
    # t_{0.975,7}, two-sided interval for the eight independent draw means.
    t_critical = 2.3646242515927844
    for metric in METRICS:
        values = {name: ({k: a[before][metric] for k, a in cells.items()},
                         {k: a[after][metric] for k, a in cells.items()})
                  for name, (before, after) in CONTRASTS.items()}
        # Compare the effects of the same packet noise with and without the learner.
        noise_effects = []
        for name in ("packet_noise_without_learner", "packet_noise_with_learner"):
            before, after = values[name]
            noise_effects.append({k: after[k] - before[k] for k in cells})
        values["packet_noise_learner_interaction"] = noise_effects
        for name, (before, after) in values.items():
            delta = {key: after[key] - before[key] for key in cells}
            draws = [mean(d for key, d in delta.items() if key[3] == seed) for seed in draw_ids]
            effect = mean(draws)
            margin = t_critical * stdev(draws) / math.sqrt(len(draws))
            summaries.append(dict(contrast=name, metric=metric, cells=len(cells), draws=len(draws),
                                  baseline=mean(before.values()), comparison=mean(after.values()),
                                  delta=effect, low=effect-margin, high=effect+margin,
                                  improved=sum(d < -1e-6 for d in delta.values()),
                                  worse=sum(d > 1e-6 for d in delta.values()), worst=max(delta.values())))
    return summaries


def check_panel(folder, stem, draw_ids, conditions):
    with (folder / f"{stem}.csv").open() as stream:
        expected = summarize(csv.DictReader(stream), draw_ids, conditions)
    with (folder / f"{stem}-paired.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    published = {(row["contrast"], row["metric"]): row for row in rows}
    if len(rows) != len(published) or len(published) != len(expected):
        raise ValueError("Incomplete or duplicate installation contrasts")
    for item in expected:
        row = published[item["contrast"], item["metric"]]
        for field in item.keys() - {"contrast", "metric"}:
            if not math.isclose(float(row[field]), item[field], abs_tol=1e-10, rel_tol=1e-10):
                raise ValueError(f"Incorrect {item['contrast']} {item['metric']} {field}")
    print(f"{stem}: {len(expected)} paired navigation summaries reproduced")


def check_evidence(folder):
    check_panel(folder, "installation-screen", DRAWS, HEIGHT_DIP)
    check_panel(folder, "installation-screen-weak", WEAK_DRAWS, WEAK_HEIGHT_DIP)
    with (folder / "stress.csv").open() as stream:
        native = {(r["case"], r["scenario"]): r for r in csv.DictReader(stream) if int(r["variant"]) == 0}
    with (folder / "installation-screen-stress.csv").open() as stream:
        stress = list(csv.DictReader(stream))
    keys = {(r["case"], r["scenario"], int(r["variant"])) for r in stress}
    if len(stress) != 28 or keys != {(*key, v) for key in native for v in (0, 11)}:
        raise ValueError("Incomplete installation stress panel")
    for row in stress:
        ref = native[row["case"], row["scenario"]]
        if any(float(row[k]) != float(ref[k]) for k in ("seed", "dip", "height")):
            raise ValueError(f"Unmatched stress fixture: {row}")
        navigation_scores(row)
        if int(row["variant"]) == 0 and any(not math.isclose(float(row[k]), float(ref[k]), abs_tol=1e-6, rel_tol=1e-6)
                                            for k in METRICS if k != "abs_signed_yaw_deg"):
            raise ValueError(f"Native stress control changed: {row}")
    print("All 14 native installation stress controls reproduced")


def main():
    from marine_mag_study import OUT, ROOT, write
    from sim_dataset import REFERENCE_CSV_SHA256_BY_CASE, digest

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4)
    panel = parser.add_mutually_exclusive_group()
    panel.add_argument("--weak", action="store_true", help="confirm the fixed threshold on independent low-wave draws at two dips")
    panel.add_argument("--stress", action="store_true", help="repeat the seven earlier stress fixtures in each direction")
    args = parser.parse_args()
    if args.stress:
        from pathlib import Path
        from tempfile import TemporaryDirectory
        from marine_mag_yaw import stress_inputs
        with TemporaryDirectory() as temporary:
            tasks = [(case, 1.5, 66.5, 3037, scenario, variant, path)
                     for case, scenario, path in stress_inputs(Path(temporary)) for variant in (0, 11)]
            with ThreadPoolExecutor(args.jobs) as pool:
                stress = list(pool.map(simulate, tasks))
        write(OUT / "installation-screen-stress.csv", stress)
        return
    stem, draws, conditions = (("installation-screen-weak", WEAK_DRAWS, WEAK_HEIGHT_DIP) if args.weak
                               else ("installation-screen", DRAWS, HEIGHT_DIP))
    tasks = []
    for case in ("following", "head"):
        for height, dip in conditions:
            folder = ROOT / ".sim-data/cases" / case
            path, = folder.glob(f"wave_data_pmstokes_H{height:.3f}_*.csv")
            if digest(path) != REFERENCE_CSV_SHA256_BY_CASE[case][path.name]:
                raise ValueError(f"Unexpected release input {path}")
            tasks.extend((case, height, dip, seed, scenario, variant, path)
                         for seed in draws for scenario in SCENARIOS for variant in VARIANTS)
    with ThreadPoolExecutor(args.jobs) as pool:
        rows = []
        for row in pool.map(simulate, tasks):
            rows.append(row)
            if len(rows) % (len(draws) * len(SCENARIOS) * len(VARIANTS)) == 0:
                print(f"{row['case']} {row['height']:g} m, dip {row['dip']:g} complete ({len(rows)}/{len(tasks)} replays)", flush=True)
    write(OUT / f"{stem}.csv", rows)
    write(OUT / f"{stem}-paired.csv", summarize(rows, draws, conditions))


if __name__ == "__main__":
    main()
