#!/usr/bin/env python3
"""Paired full-rate TFG replays. Keep transient logs outside the evidence tree."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from sim_dataset import ROOT, digest, REFERENCE_CSV_SHA256_BY_CASE
from tfg_shared_modes_degradation import (DIAGNOSTIC_HIGH_HEIGHT_M, DIAGNOSTIC_SOURCE_HEIGHT_M,
    ROUNDTRIP_DURATION_SEC, ROUNDTRIP_WAVE_NAME, make_roundtrip_wave, ov)

BINARY = Path(os.environ.get('MARINE_MAG_BINARY', ROOT/'tests/kalman_tfg_shared_modes/magnetic_study/build/magnetic-study'))
OUT = ROOT / 'reports/results/marine_mag'
CAL = 't bx by bz px py pz applied_matrix_err matrix_err field_err dip_err accepted missed rejected first_applied info_hard info_scale scale_open bias_cov_trace bias_nees min_eig'.split()
NAV = 'yaw_rms_deg roll_rms_deg pitch_rms_deg disp_x_rms_m disp_y_rms_m disp_z_rms_m disp_3d_rms_m gyro_bias_3d_rms_radps accel_bias_3d_rms_mps2'.split()
ARMS = tuple(range(-1,9))


def plan(stage):
    """A split is a sensor/installation draw; seas within a draw are paired."""
    for case in ('following','head'):
        if stage=='diagnosis':
            for scenario in ('known','ideal','clean','hard','diagonal','cross','alignment','combined',
                             'dip_error','tilt_kick','yaw_kick','uncoupled','latency','burst','interval'):
                for v in (-1,0):
                    yield case,'jonswap','1.500',0,scenario,66.5,v
        if stage=='development':
            for seed in (13,29):
                for height,s in [('0.270','clean'),('1.500','hard'),('4.000','combined')]:
                    for v in (*ARMS,9):
                        yield case,'jonswap',height,seed,s,66.5,v
        if stage=='validation':
            for seed in (101,307,1009,2027):
                # Fixed blocked challenge panel, with previously unseen spectra,
                # inclinations, installation draws and corruption/noise draws.
                for height,s,dip in [('0.270','clean',45),('1.500','combined',66.5),
                                     ('8.500','combined',76),('4.000','burst',55)]:
                    for v in ARMS:
                        yield case,'pmstokes',height,seed,s,dip,v
        if stage=='stress':
            for s in ('engine','mount','outlier','latency','roundtrip','weak','turn'):
                for v in ARMS:
                    yield case,'jonswap','1.500',3037,s,66.5,v


def write(path,rows):
    if not rows:
        return
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def special(data,case,cache,scenario):
    """Rebuild all dependent channels using project helpers and the release clock."""
    folder=cache/'inputs'/case/scenario
    folder.mkdir(parents=True,exist_ok=True)
    low=ov.find_default_input(data,'1.500','50.710')
    path=folder/(ROUNDTRIP_WAVE_NAME if scenario=='roundtrip' else low.name)
    if digest(low) != REFERENCE_CSV_SHA256_BY_CASE[case][low.name]:
        raise ValueError(f'Unexpected source for regenerated maneuver: {low}')
    # underway_reference selected_case() resolves OCEAN_IMU_SIM_CASE at call time
    os.environ['OCEAN_IMU_SIM_CASE']=case
    columns,values=ov.read_wave_csv(low,ROUNDTRIP_DURATION_SEC)
    if scenario=='roundtrip':
        hi=ov.find_default_input(data,'8.500','202.839')
        if digest(hi) != REFERENCE_CSV_SHA256_BY_CASE[case][hi.name]:
            raise ValueError(f'Unexpected high-sea source: {hi}')
        _,high=ov.read_wave_csv(hi,ROUNDTRIP_DURATION_SEC)
        values=make_roundtrip_wave(ov,columns,values,high,11,
                                  DIAGNOSTIC_HIGH_HEIGHT_M/DIAGNOSTIC_SOURCE_HEIGHT_M)
    else:
        from scipy.integrate import cumulative_trapezoid
        from underway_reference import track,SPEED_MPS
        t=values[:,columns.index('time')]
        p,v,a,heading=track(t,case)
        mean=0 if case=='following' else 180
        yaw=mean+np.sin(2*np.pi*t/20)*(1 if scenario=='weak' else 3)
        if scenario=='turn':
            yaw+=30*(np.tanh((t-550)/12)-np.tanh((t-630)/12))/2
        angle=np.deg2rad(yaw)
        vel=SPEED_MPS*np.column_stack((np.cos(angle),np.sin(angle),np.zeros_like(t)))
        acc=np.gradient(vel,t,axis=0,edge_order=2)
        pos=cumulative_trapezoid(vel,t,axis=0,initial=0)
        for prefix,old,new in [('disp',p,pos),('vel',v,vel),('acc',a,acc)]:
            idx=[columns.index(prefix+'_'+c) for c in 'xyz']
            values[:,idx]+=new-old
        values[:,columns.index('yaw_deg')]+=yaw-heading
        values=ov.rebuild_body_imu(columns,values)
        # Score-only correction to W3D's fixed release-track subtraction.
        # Keep p, v, acceleration, attitude and IMU physically consistent in CSV.
    with tempfile.TemporaryDirectory(dir=folder) as tmp:
        staged=Path(tmp)/path.name
        ov.write_wave_csv(staged,columns,values)
        if scenario!='roundtrip':
            np.savetxt(staged.with_suffix('.track'),(pos-p)[:,:2],fmt='%.12g')
            staged.with_suffix('.track').replace(path.with_suffix('.track'))
        staged.replace(path)
    return path


def physical_scenario(scenario, env):
    """Apply the simulator fixture for a named physical stress case."""
    physical=scenario
    if scenario in ('engine','mount','interval','roundtrip'):
        physical='combined'
    if scenario=='interval': env['W3D_CSM_GYRO_INPUT']='interval'
    if scenario=='mount': env.update(W3D_IMU_LEVER_ARM_M='1,0,0',W3D_IMU_LEVER_ARM_MODEL='none')
    if scenario=='engine': env['W3D_ENGINE_RPM']='1600'
    return physical


def parse(stdout, allow_intervention_failure=False):
    metrics={}
    for line in stdout.splitlines():
        if line.startswith(('VALIDATION_METRICS ','VALIDATION_METRICS_MEANS ')):
            fields=dict(x.split('=',1) for x in line.split()[1:])
            if fields['segment']=='window':
                metrics.update(fields)
    if int(metrics['samples'])!=180000:
        raise ValueError('incomplete final 900 seconds')
    nav={k:float(metrics[k]) for k in NAV+['yaw_mean_deg']}
    nav['centered_deg']=np.sqrt(max(0,nav['yaw_rms_deg']**2-nav['yaw_mean_deg']**2))
    def rows(tag):
        return np.array([[float(v) for v in line.split()[1:]] for line in stdout.splitlines() if line.startswith(tag+' ')])
    c=rows('CAL')
    y=rows('YAW')
    if c.shape!=(239,len(CAL)) or y.shape!=(240,4) or not np.all(np.isfinite(c)) or not np.all(np.isfinite(y)):
        raise ValueError('nonfinite or missing diagnostics')
    fields=dict(zip(CAL,c.T))
    for k in ['bx','by','bz','px','py','pz','applied_matrix_err','matrix_err','field_err','dip_err']:
        nav[k+'_rms']=np.sqrt(np.mean(fields[k][60:]**2))
    for k in ['accepted','missed','rejected','first_applied','info_hard','info_scale','scale_open','bias_cov_trace']:
        nav[k]=fields[k][-1]
    nav['bias_nees_late']=np.mean(fields['bias_nees'][60:])
    nav['min_eig']=min(fields['min_eig'])
    nav['max_live_step_deg']=max(y[60:,3])
    for label,start,end in [('pre',540,600),('event',600,610),('recovery',610,640),('late',1140,1200),
                            ('turn',550,630),('after_turn',630,710)]:
        nav[label+'_yaw_deg']=np.sqrt(np.mean(y[(y[:,0]>start)&(y[:,0]<=end),2]**2))
    failed=[k for k,v in nav.items() if not np.isfinite(v)]
    if failed and not allow_intervention_failure:
        raise ValueError('nonfinite metrics: '+str(failed))
    for k in failed:
        nav[k]=None
    nav['failed_metrics']=','.join(failed)
    return nav,np.column_stack((y[:len(c)],c[:,1:])).tolist()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['diagnosis','development','validation','stress'])
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--out',type=Path,default=OUT)
    args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    identity=digest(BINARY)
    scorer=digest(Path(__file__))
    conditions=list(plan(args.stage))
    records={}
    for case,spectrum,height,_seed,scenario,_dip,_variant in conditions:
        key=case,spectrum,height,scenario
        if key in records:
            continue
        data=ROOT/'.sim-data/cases'/case
        path=(special(data,case,args.cache,scenario) if scenario in ('turn','weak','roundtrip')
              else next(data.glob(f'wave_data_{spectrum}_H{height}_*.csv')))
        sha=digest(path)
        if scenario not in ('turn','weak','roundtrip') and sha!=REFERENCE_CSV_SHA256_BY_CASE[case][path.name]:
            raise ValueError(f'Unexpected release input {path}')
        track_sha=digest(path.with_suffix('.track')) if scenario in ('turn','weak') else None
        records[key]=path,sha,track_sha
    def run(condition):
        case,spectrum,height,seed,scenario,dip,v=condition
        path,sha,track_sha=records[case,spectrum,height,scenario]
        ck=hashlib.sha256(json.dumps([identity,sha,track_sha,scorer,condition]).encode()).hexdigest()
        checkpoint=args.cache/(ck+'.json')
        key=dict(case=case,spectrum=spectrum,height=height,seed=seed,scenario=scenario,dip=dip,variant=v)
        if not checkpoint.exists():
            env={k:val for k,val in os.environ.items() if not k.startswith(('W3D_','OCEAN_IMU_SIM_CASE'))}
            env['OCEAN_IMU_SIM_CASE']=case
            if seed: env['W3D_SEED']=str(seed)
            physical = physical_scenario(scenario, env)
            run=subprocess.run([BINARY,path,str(v),physical,str(dip)],env=env,capture_output=True,text=True)
            if run.returncode: raise RuntimeError(str(condition)+' '+run.stderr[-1000:]+run.stdout[-800:])
            try:
                scores,trace=parse(run.stdout,scenario=="uncoupled")
            except ValueError as e:
                (args.cache/(ck+'.failed.log')).write_text(run.stdout)
                raise ValueError(f'{condition}: {e}') from e
            with tempfile.NamedTemporaryFile(mode='w',dir=args.cache,delete=False) as stream:
                json.dump([scores,trace],stream,allow_nan=False)
            Path(stream.name).replace(checkpoint)
        scores,trace=json.loads(checkpoint.read_text())
        scores.setdefault('failed_metrics','')
        # Truth-qualified convergence is a scoring outcome only. Clean controls
        # may qualify immediately even with an identity installation correction.
        history=np.asarray(trace)
        def column(k):
            return history[:,4+CAL[1:].index(k)]
        hard_error=np.sqrt(sum(column(k)**2 for k in ('bx','by','bz')))
        usable=(history[:,0]>=180)&(hard_error<=1.)&(column('applied_matrix_err')<=.03)
        held=np.flatnonzero(np.convolve(usable.astype(int),np.ones(13,dtype=int),'valid')==13)
        scores['usable_s']=int(history[held[0],0]+60) if len(held) else -1
        scores['bias_cover_95']=np.mean(column('bias_nees')[60:]<=7.814728) if v!=-1 else None
        print(args.stage,case,height,seed,scenario,v,round(scores['yaw_rms_deg'],4),flush=True)
        return dict(**key,**scores),[dict(**key,**dict(zip(['t','yaw_mean','yaw_rms','peak_step']+CAL[1:],row))) for row in trace]
    with ThreadPoolExecutor(args.jobs) as pool: results=list(pool.map(run,conditions))
    write(args.out/(args.stage+'.csv'),[r for r,_ in results])
    # Uniform five-second diagnostics on one matched family; full logs are reproducible.
    write(args.out/(args.stage+'-trace.csv'),[q for r,t in results if
        r['case']=='following' and ((args.stage=='diagnosis' and r['scenario'] in ('hard','known','combined') and r['variant']==0)
        or (args.stage=='development' and r['height']=='4.000' and r['seed']==13 and r['variant'] in (0,2,6,7))
        or (args.stage=='stress' and r['scenario'] in ('burst','turn','roundtrip') and r['variant'] in (0,6))
        or (args.stage=='validation' and r['scenario']=='burst' and r['seed']==101 and r['variant'] in (0,6,7)))
        for q in t[::2]])
    write(args.out/(args.stage+'-inputs.csv'),[dict(case=k[0],spectrum=k[1],height=k[2],scenario=k[3],file=p.name,sha256=h,
          **({'track_sha256':track_sha or ''} if args.stage=='stress' else {}))
          for k,(p,h,track_sha) in sorted(records.items())])


if __name__=='__main__':
    main()
