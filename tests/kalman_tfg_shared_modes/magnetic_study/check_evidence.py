"""Fail closed on missing, unpaired, nonfinite or numerically inconsistent holdout rows."""
import csv
import math
import sys
from pathlib import Path
from statistics import mean

root=Path(__file__).resolve().parents[3]
with (root/'reports/results/marine_mag/validation.csv').open() as f:
    rows=list(csv.DictReader(f))
groups={}
for row in rows:
    k=tuple(row[n] for n in ('case','height','seed','scenario','dip'))
    variant=int(row['variant'])
    assert variant not in groups.setdefault(k,{})
    groups[k][variant]=row
    for field in ('yaw_rms_deg','yaw_mean_deg','centered_deg','disp_3d_rms_m','min_eig'):
        assert math.isfinite(float(row[field])), (k,variant,field)
    assert float(row['min_eig'])>=-1e-5
    assert math.isclose(float(row['yaw_rms_deg'])**2,
        float(row['yaw_mean_deg'])**2+float(row['centered_deg'])**2,abs_tol=1e-5)
assert len(groups)==32
for key,arms in groups.items():
    assert set(arms)==set(range(-1,9)),key
    for v in (1,5):
        for field in ('yaw_rms_deg','yaw_mean_deg','roll_rms_deg','pitch_rms_deg','disp_3d_rms_m','bx_rms'):
            assert abs(float(arms[0][field])-float(arms[v][field]))<=1e-7,(key,v,field)
    assert abs(float(arms[0]['matrix_err_rms'])-float(arms[7]['matrix_err_rms']))<1e-7
with (root/'reports/results/marine_mag/paired.csv').open() as f:
    yaw_rows=[r for r in csv.DictReader(f) if r['metric']=='yaw_rms_deg']
    expected=set(range(-1,9))
    assert len(yaw_rows)==len(expected) and {int(r['variant']) for r in yaw_rows}==expected
    scores={int(r['variant']):float(r['delta']) for r in yaw_rows}
for v,delta in scores.items():
    actual=mean(float(g[v]['yaw_rms_deg'])-float(g[0]['yaw_rms_deg']) for g in groups.values())
    assert abs(actual-delta)<1e-6
print(f'{len(groups)} complete paired conditions, {len(rows)} finite scores; covariance and control identities pass')


with (root / "reports/results/marine_mag/stress-inputs.csv").open() as stream:
    fixtures = list(csv.DictReader(stream))
tracks = [r for r in fixtures if r["scenario"] in ("turn", "weak")]
assert len(tracks) == 4
for row in tracks:
    assert len(bytes.fromhex(row["track_sha256"])) == 32

with (root / "reports/results/marine_mag/quarantine-only.csv").open() as stream:
    quarantine_rows = list(csv.DictReader(stream))
assert len(quarantine_rows) == len(groups)
seen = set()
for row in quarantine_rows:
    key = tuple(row[name] for name in ("case", "height", "seed", "scenario", "dip"))
    assert key not in seen and key in groups
    seen.add(key)
    native = groups[key][0]
    assert float(row["matrix_err_rms"]) == float(native["matrix_err_rms"])
    assert math.isfinite(float(row["yaw_rms_deg"])), (key, "nonfinite quarantine yaw")
    if row["scenario"] != "burst":
        assert float(row["rejected"]) == 0 and float(row["yaw_rms_deg"]) == float(native["yaw_rms_deg"])


# Recompute all published holdout effect means for the two sustained-yaw studies.
for study, stage, factor, expected_seeds in (
    ("yaw-weight", "validation", "sigma_factor", {421,829,1609,3517}),
    ("gyro-bias", "validation", "gyro_bias_factor", {547,1151,2371,4099}),
    ("gyro-bias", "replication", "gyro_bias_factor", {5003,5011,5021,5023,5039,5051,5059,5077,5081,5087,5099,5101}),
):
    folder = root / "reports/results/marine_mag"
    with (folder / f"{study}-{stage}.csv").open() as stream:
        arms = {}
        for row in csv.DictReader(stream):
            key = tuple(row[name] for name in ("case","height","seed","scenario","dip"))
            scale = float(row[factor])
            assert scale not in arms.setdefault(key, {})
            arms[key][scale] = row
            assert float(row["min_eig"]) >= -1e-5 and not row["failed_metrics"]
    for seed in expected_seeds:
        assert {(k[0],float(k[1])) for k in arms if int(k[2])==seed} == {
            (case,height) for case in ("following","head") for height in (.27,1.5,4.,8.5)}
    assert len(arms)==8*len(expected_seeds)
    scales = set(next(iter(arms.values())))
    assert len(scales) == 2 and 1.0 in scales and all(set(cell)==scales for cell in arms.values())
    stats_file = f"{study}-paired.csv" if stage=="validation" else f"{study}-{stage}-paired.csv"
    with (folder / stats_file).open() as stream:
        published = [r for r in csv.DictReader(stream) if r.get("stage",stage) == stage]
    metrics = {"yaw_rms_deg","abs_signed_yaw_deg","centered_deg","roll_rms_deg","pitch_rms_deg","disp_3d_rms_m"}
    assert len(published)==len(scales)*len(metrics)
    assert {(float(r[factor]),r["metric"]) for r in published}=={(s,m) for s in scales for m in metrics}
    for row in published:
        scale, metric = float(row[factor]), row["metric"]
        delta = mean(float(cell[scale][metric])-float(cell[1.0][metric]) for cell in arms.values())
        assert math.isclose(float(row["delta"]),delta,abs_tol=1e-8)
    print(f"{study} {stage}: {len(arms)} complete paired cells and all published mean differences reproduced")


sys.path.insert(0, str(root / "tools"))
from marine_mag_attribution import check_evidence

check_evidence(root / "reports/results/marine_mag")
