"""Paired, record-balanced analysis of the RAO coefficient response surface."""
from __future__ import annotations
import argparse
import csv
import json
import hashlib
import math
import re
from itertools import product
from collections import defaultdict
from pathlib import Path
import numpy as np

METRICS = ('disp_z_pct_hs','disp_3d_rms_m','roll_rms_deg','pitch_rms_deg',
           'yaw_rms_deg','accel_bias_3d_rms_mps2','gyro_bias_3d_rms_radps')
WEIGHTS = np.array([2.,2.,.5,.5,.25,.75,.5])
KEYS = ('TFG_TAU_COEFF','TFG_SIGMA_COEFF','TFG_R_S_X_FACTOR','TFG_R_S_Y_FACTOR')

def parameters(env):
    return [float(env[k]) for k in KEYS]

def analyse(rows: list[dict]) -> dict:
    grouped = defaultdict(dict)
    for row in rows:
        key = (row['input'],str(row['seed']))
        if key in grouped[row['config']]:
            raise ValueError(f'duplicate case {row["config"]} {key}')
        grouped[row['config']][key] = row
    base = grouped['baseline']
    cases = sorted(base)
    if len({key[0] for key in cases}) != 8:
        raise ValueError('baseline does not contain all eight RAO records')
    inputs = {k[0] for k in cases}
    seeds = {k[1] for k in cases}
    if set(cases) != set(product(inputs, seeds)):
        raise ValueError('incomplete baseline record/seed cross product')
    if any(base[k].get('error') is not None or
           not isinstance(base[k]['metrics'].get(m),(int,float)) or
           not math.isfinite(base[k]['metrics'][m]) or base[k]['metrics'][m] < 0
           for k in cases for m in METRICS):
        raise ValueError('baseline contains invalid physical metrics')
    output = []
    for name, data in grouped.items():
        report = dict(config=name, parameters=parameters(next(iter(data.values()))['env']))
        missing = set(cases)-set(data)
        unexpected = set(data)-set(cases)
        if unexpected:
            raise ValueError(f'{name}: unexpected record/seed cases: {sorted(unexpected)}')
        invalid = [key for key in cases if key in data and
                   (data[key].get('error') is not None or any(not isinstance(data[key].get('metrics',{}).get(m),(int,float)) or
                       not math.isfinite(data[key]['metrics'][m]) or
                       data[key]['metrics'][m]<0 for m in METRICS))]
        report.update(missing_cases=sorted(missing), invalid_cases=invalid)
        if missing or invalid:
            report.update(score=None, admissible=False)
            output.append(report)
            continue
        current = np.array([[data[k]['metrics'][m] for m in METRICS] for k in cases])
        reference = np.array([[base[k]['metrics'][m] for m in METRICS] for k in cases])
        ratios = np.maximum(current,1e-15)/np.maximum(reference,1e-15)
        logratios = np.log(ratios)
        new_failures = [k for k in cases if data[k]['exit_code']!=0 and base[k]['exit_code']==0]
        extra_violations = sum(max(0,len(data[k]['violations'])-len(base[k]['violations'])) for k in cases)
        raw_score = float(np.mean(logratios @ WEIGHTS/WEIGHTS.sum()))
        penalty = .025*extra_violations
        changes = {}
        for j,m in enumerate(METRICS):
            record_changes = {}
            record_mean_changes = {}
            seed_changes = {}
            for input_name in sorted({k[0] for k in cases}):
                idx=[i for i,k in enumerate(cases) if k[0]==input_name]
                record_mean_changes[input_name]=float(100*(current[idx,j].mean()/reference[idx,j].mean()-1))
                record_changes[input_name]=float(100*(np.sqrt(np.mean(current[idx,j]**2))/np.sqrt(np.mean(reference[idx,j]**2))-1))
            for seed in sorted({k[1] for k in cases}):
                idx=[i for i,k in enumerate(cases) if k[1]==seed]
                seed_changes[seed]=float(100*(np.sqrt(np.mean(current[idx,j]**2))/np.sqrt(np.mean(reference[idx,j]**2))-1))
            changes[m]=dict(
                baseline_mean=float(reference[:,j].mean()), candidate_mean=float(current[:,j].mean()),
                mean_rms_change_pct=float(100*(current[:,j].mean()/reference[:,j].mean()-1)),
                pooled_change_pct=float(100*(np.sqrt(np.mean(current[:,j]**2))/np.sqrt(np.mean(reference[:,j]**2))-1)),
                candidate_pooled_rms=float(np.sqrt(np.mean(current[:,j]**2))),
                baseline_pooled_rms=float(np.sqrt(np.mean(reference[:,j]**2))),
                balanced_geomean_change_pct=float(100*np.expm1(logratios[:,j].mean())),
                paired_wins=int((ratios[:,j]<1-1e-9).sum()), paired_cases=len(cases),
                record_wins=sum(v<-1e-7 for v in record_changes.values()), records=8,
                per_record_change_pct=record_changes, per_record_mean_rms_change_pct=record_mean_changes,
                mean_rms_record_wins=sum(v<-1e-7 for v in record_mean_changes.values()), per_seed_change_pct=seed_changes,
                worst_paired_change_pct=float(100*(ratios[:,j].max()-1)))
        regime_scores={}
        for label,lo,hi in [('low',0,1.5),('high',4,9)]:
            idx=[i for i,k in enumerate(cases) if lo<=float(re.search(r'_H([\d.]+)_',k[0])[1])<=hi]
            regime_scores[label]=float(np.mean(logratios[idx] @ WEIGHTS/WEIGHTS.sum()))
        report.update(score=raw_score+penalty,raw_score=raw_score,penalty=penalty,
                      admissible=True,new_gate_failure_cases=new_failures,extra_gate_violations=extra_violations,
                      failed_cases=sum(data[k]['exit_code']!=0 for k in cases),
                      baseline_failed_cases=sum(base[k]['exit_code']!=0 for k in cases),
                      metrics=changes,regime_scores=regime_scores)
        output.append(report)
    output.sort(key=lambda x:(x['score'] is None,x['score'] if x['score'] is not None else 0,x['config']))
    valid=[r for r in output if r['admissible']]
    pareto=[]
    for r in valid:
        vector=np.array([r['metrics'][m]['balanced_geomean_change_pct'] for m in METRICS])
        if not any(np.all((other:=np.array([q['metrics'][m]['balanced_geomean_change_pct'] for m in METRICS]))<=vector) and
                   np.any(other<vector) for q in valid if q is not r):
            pareto.append(r['config'])
    points=np.array([r['parameters'] for r in valid])
    bounds=dict(zip(KEYS,zip(points.min(axis=0).tolist(),points.max(axis=0).tolist())))
    boundary={}
    for r in valid[:10]:
        boundary[r['config']]=[]
        for j,x in enumerate(r['parameters']):
            for side,b in [('low',points[:,j].min()),('high',points[:,j].max())]:
                if math.isclose(x,b,rel_tol=1e-9):boundary[r['config']].append(KEYS[j]+':'+side)
    joint = [r for r in valid if r['config'].startswith('joint_')]
    joint_contacts = {}
    if joint:
        joint_points = np.array([r['parameters'] for r in joint])
        for r in joint[:5]:
            joint_contacts[r['config']] = [
                KEYS[j] + ':' + side
                for j, value in enumerate(r['parameters'])
                for side, bound in [('low', joint_points[:, j].min()),
                                    ('high', joint_points[:, j].max())]
                if math.isclose(value, bound, rel_tol=1e-9)]
    features=np.column_stack((np.log(points[:,0]),np.log(points[:,1]),
                             .5*np.log(points[:,2]*points[:,3]),np.log(points[:,2]/points[:,3])))
    active=[i for i in range(4) if np.ptp(features[:,i])>1e-10]
    x=features[:,active]
    mu=x.mean(axis=0); scale=x.std(axis=0); x=(x-mu)/scale
    columns=[np.ones(len(x))]+[x[:,i] for i in range(x.shape[1])]
    labels=['intercept']+[f'linear_{i}' for i in active]
    for i in range(x.shape[1]):
        for j in range(i,x.shape[1]):
            columns.append(x[:,i]*x[:,j]);labels.append(f'interaction_{active[i]}_{active[j]}')
    design=np.column_stack(columns)
    responses=np.array([r['raw_score'] for r in valid])
    metric_surfaces={}
    for metric in METRICS:
        response=np.array([math.log1p(r['metrics'][metric]['balanced_geomean_change_pct']/100) for r in valid])
        coeff,_,_,_=np.linalg.lstsq(design,response,rcond=None)
        residual=design@coeff-response
        metric_surfaces[metric]={'coefficients':dict(zip(labels,coeff.tolist())),
            'r_squared':1-float(residual@residual)/max(float(np.sum((response-response.mean())**2)),1e-15)}
    coefficients,_,rank,singular=np.linalg.lstsq(design,responses,rcond=None)
    predicted=design@coefficients
    r2=1-float(np.sum((predicted-responses)**2))/max(float(np.sum((responses-responses.mean())**2)),1e-15)
    return dict(analysis_producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                aggregation={'pooled_change_pct':'sqrt(mean(case_RMS**2)) ratio',
                    'mean_rms_change_pct':'mean(case_RMS) ratio',
                    'balanced_geomean_change_pct':'equal-case geometric mean of paired ratios'},
                ranking=output,pareto=pareto,bounds=bounds,top_boundary_contacts=boundary,coarse_joint_face_contacts=joint_contacts,
                regime_winners={g:min(valid,key=lambda r:r['regime_scores'][g])['config'] for g in ('low','high')},
                surface_diagnostic=dict(feature_names=['log_tau','log_sigma','log_geometric_xy','log_x_over_y'],
                    active_features=active,normalization_mean=mu.tolist(),normalization_scale=scale.tolist(),
                    coefficients=dict(zip(labels,coefficients.tolist())),rank=int(rank),r_squared=r2,
                    condition_number=float(singular[0]/singular[-1]),metric_surfaces=metric_surfaces,
                    caveat='Exploratory quadratic summary, not a substitute for simulations or a global optimality certificate.'))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('runs',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();result=analyse(json.loads(args.runs.read_text()));args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    with (args.output/'response_surface.csv').open('w',newline='') as handle:
        names=['config',*KEYS,'score','raw_score','extra_gate_violations',*METRICS]
        writer=csv.DictWriter(handle,fieldnames=names);writer.writeheader()
        for row in result['ranking']:
            flat=dict(config=row['config'],**dict(zip(KEYS,row['parameters'])),score=row.get('score'),
                      raw_score=row.get('raw_score'),extra_gate_violations=row.get('extra_gate_violations'))
            flat.update({m:row.get('metrics',{}).get(m,{}).get('balanced_geomean_change_pct') for m in METRICS})
            writer.writerow(flat)
    print(json.dumps({k:v for k,v in result.items() if k!='ranking'},indent=2))
    print('Top candidates:',[(r['config'],r['score']) for r in result['ranking'][:10]])
if __name__=='__main__':main()
