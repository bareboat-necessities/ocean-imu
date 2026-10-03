# ruff: noqa: F401, F811
"""Root parameterization for conditional OU-III source-uniform history covers."""
from __future__ import annotations
import json
from pathlib import Path
from .reachable_history_enclosure import HistoryCell,HistoryInterval
ROOT=Path(__file__).resolve().parents[3]
C=json.loads((ROOT/"tools/stability/ou3_theorem/constants.json").read_text())

def root_cells(horizon_s=60.):
 """Finite-dimensional causal generator coordinates, not generated filter outputs.

 Piecewise history knots are proof coordinates. Temporal constraints in
 temporal_source_domain reject inadmissible combinations; they are not reset at
 word boundaries.
 """
 if horizon_s not in (60.,100.):raise ValueError("certified horizons 60/100 s")
 imu=C["imu_bias"];m=C["marine_motion"]
 coords=[]
 knots=[0.,30.,60.] if horizon_s==60 else [0.,30.,60.,90.,100.]
 for t in knots:
  for a in range(3):
   coords += [(f"slow_accel_{a}_{t}",HistoryInterval(-imu["B_a_s_mps2"],imu["B_a_s_mps2"])),
              (f"slow_gyro_{a}_{t}",HistoryInterval(-imu["B_g_s_rad_s"],imu["B_g_s_rad_s"])),
              (f"fast_accel_primitive_{a}_{t}",HistoryInterval(-imu["fast_accel_accumulation_cap_mps"],imu["fast_accel_accumulation_cap_mps"])),
              (f"fast_gyro_primitive_{a}_{t}",HistoryInterval(-imu["fast_gyro_accumulation_cap_rad"],imu["fast_gyro_accumulation_cap_rad"])),
              (f"physical_p_{a}_{t}",HistoryInterval(-m["P_max_m"],m["P_max_m"])),
              (f"physical_v_{a}_{t}",HistoryInterval(-m["V_max_mps"],m["V_max_mps"]))]
 # Attitude/magnetic histories are parameterized by delivered vector knots,
 # not Mahony/reference states. Causal propagation generates those states.
 for t in knots:
  for a in range(3):
   coords += [(f"gravity_dir_{a}_{t}",HistoryInterval(-1,1)),
              (f"mag_field_{a}_{t}",HistoryInterval(-C["magnetic_service"]["field_norm_max_uT"],C["magnetic_service"]["field_norm_max_uT"]))]
 return (HistoryCell(tuple(coords),f"source-uniform-{int(horizon_s)}s-root"),)
