"""Read-only 17-s held-BA LIN variational diagnostic.

Records the literal shipping correction/prediction factors during a held-BA
window ending before first BA activation. AG tangent columns are treated as
inputs; this diagnostic propagates only the 12 LIN root columns.
"""
from pathlib import Path
import json, subprocess, tempfile
import numpy as np
from .construction_history_diagnostic import REPO

NL=12
OFF_L=6

def driver_source():
    s=(REPO/"tools/stability/ag_readout_source.cpp").read_text()
    # Record from Live onward so Python can select an exact held 17-s suffix.
    s=s.replace("if (k==45001) {\n            if (active<0 || k-active<3400) return 2;\n            root=matrix_json(filter.raw().mekf().covariance_full());\n            recording=true;\n        }",
                "if (filter.isLive()) recording=true;")
    s=s.replace(" || applied!=8","")
    return s

def lin_factor(ev):
    if ev["kind"]=="prediction":
        return np.asarray(ev["F_LIN"],float)
    if ev["kind"]=="correction":
        H=np.asarray(ev["H"],float)
        K=np.asarray(ev["K"],float)
        # Full H is 3 x NX; LIN block is columns 6:18.
        HL=H[:,OFF_L:OFF_L+NL]
        KL=K[OFF_L:OFF_L+NL,:]
        return np.eye(NL)-KL@HL
    # covariance sync and attitude reset do not change LIN mean directly.
    return np.eye(NL)

def run(eigen, mode="wave"):
    with tempfile.TemporaryDirectory(prefix="ou3-held-lin-") as d:
        src=Path(d)/"r.cpp"; exe=Path(d)/"r"; src.write_text(driver_source())
        subprocess.run(["g++","-O2","-std=c++20","-I"+str(REPO/"src"),"-isystem",str(eigen),str(src),"-o",str(exe)],check=True)
        native=json.loads(subprocess.check_output([str(exe),mode],text=True))
    live=native["live_step"]; active=native["active_step"]
    if active<0 or live<0 or active-live<3400:
        raise RuntimeError(f"no 17-s held window: live={live} active={active}")
    # Events are chronological from Live. Select final 17 s before BA activation.
    # Each prediction carries one sample; use exactly last 3400 predictions and all
    # interleaved operations between their first and the activation boundary.
    events=native["events"]
    pred_idx=[i for i,e in enumerate(events) if e["kind"]=="prediction"]
    cut=pred_idx[-3400]
    chosen=events[cut:]
    M=np.eye(NL)
    for ev in chosen:
        M=lin_factor(ev)@M
    eigvals=np.linalg.eigvals(M)
    rho=float(np.max(np.abs(eigvals)))
    smax=float(np.linalg.svd(M,compute_uv=False)[0])
    return {"mode":mode,"live_step":live,"active_step":active,
            "held_window_predictions":3400,"held_window_s":17.0,
            "event_count":len(chosen),"rho_euclidean":rho,
            "sigma_max_euclidean":smax,
            "eigenvalues":[[float(z.real),float(z.imag)] for z in eigvals],
            "qualification":"FINITE_CARRIED_HELD_BA_DIAGNOSTIC_ONLY",
            "source_uniform_verified":False}

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(); p.add_argument("--eigen",type=Path,default=Path("/usr/include/eigen3")); p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    out={"quiet":run(a.eigen,"0"),"wave":run(a.eigen,"wave")}
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2,sort_keys=True))
