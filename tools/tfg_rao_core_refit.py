#!/usr/bin/env python3
"""Focused TFG RAO refit with independent selection and sealed holdout phases.

The explicit extension/refinement plan follows the completed 2,088-case broad
screen. No completed generic or attitude/bias campaign is repeated. The report
interface and candidate freezing are retained; all replays are resumable.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = json.loads((ROOT/'tools/tfg_rao_refit_protocol.json').read_text())
KEYS = tuple(PROTOCOL['coefficient_keys'])
W = PROTOCOL['physical_metric_weights']
METRICS = tuple(W)
LEGACY = (1., .8, 1.15, 1.15)
BASEENV = {'TFG_TUNING':'adaptive', 'TFG_AW_COV_SYNC':'1',
           'SF_GYRO_BIAS_RW_VAR':'3.75e-10'}
SEEDS = {'screen':PROTOCOL['training_seeds'], 'extend':PROTOCOL['training_seeds'],
         'split':PROTOCOL['training_seeds'],
         'refine':PROTOCOL['refinement_seeds'], 'holdout':PROTOCOL['sealed_holdout_seeds']}


def arm(name, tau=1., sigma=.8, x=1.15, y=1.15):
    return {'name':name, 'env':dict(BASEENV, **dict(zip(KEYS,map(str,(tau,sigma,x,y)))))}


def coeffs(env):
    return tuple(float(env.get(k,v)) for k,v in zip(KEYS,LEGACY))


def cases(rows):
    grouped = {}
    for row in rows:
        key = (row['input'],str(row['seed']))
        group = grouped.setdefault(row['config'],{})
        if key in group:
            raise ValueError(f'duplicate case: {row["config"]} {key}')
        group[key] = row
    return grouped


def worst(row):
    return max([v['ratio'] for v in row['violations'] if v['ratio'] is not None],default=1.)


def rank(rows):
    """Predeclared paired score; usable also for report record/seed subsets."""
    grouped = cases(rows)
    base = grouped['baseline']
    result = []
    for name,group in grouped.items():
        good = set(group)==set(base)
        logs = {m:[] for m in METRICS}
        extra = 0
        for key,row in group.items():
            reference = base.get(key)
            if reference is None:
                good = False
                continue
            # Net totals can conceal a new violation on a previously good
            # history. The predeclared penalty is per paired case.
            extra += max(0,len(row['violations'])-len(reference['violations']))
            if row.get('error') or reference.get('error'):
                good = False
            for metric in METRICS:
                value = row['metrics'].get(metric)
                denominator = reference['metrics'].get(metric)
                if (not isinstance(value,(int,float)) or not isinstance(denominator,(int,float))
                        or not math.isfinite(value) or not math.isfinite(denominator)
                        or value<0 or denominator<=0):
                    good = False
                    continue
                logs[metric].append(math.log(max(value,1e-15)/denominator))
        ratios = {m:math.exp(sum(v)/len(v)) if v else None for m,v in logs.items()}
        raw = sum(W[m]*sum(v)/len(v) for m,v in logs.items())/sum(W.values()) if good else math.inf
        penalty = .025*extra
        result.append(dict(config=name,coefficients=dict(zip(KEYS,coeffs(next(iter(group.values()))['env']))),
                           eligible=good,cases=len(group),expected_cases=len(base),
                           score=raw+penalty,raw_score=raw,penalty=penalty,
                           gate_failures=sum(r['exit_code']!=0 for r in group.values()),
                           violations=sum(len(r['violations']) for r in group.values()),ratios=ratios))
    return sorted(result,key=lambda r:(r['score'],r['config']))


def screen():
    configs = [arm('baseline')]
    for v in (.45,.6,.75,.9,1.1,1.3,1.55,1.85):
        configs.append(arm(f'tau_{v}',tau=v))
    for v in (.35,.5,.65,.95,1.15,1.4):
        configs.append(arm(f'sigma_{v}',sigma=v))
    for v in (.45,.65,.85,1.,1.3,1.55,1.85,2.2):
        configs.append(arm(f'rsxy_{v}',x=v,y=v))
    for t,s,xy in itertools.product((.65,.9,1.15,1.45),(.5,.75,1.,1.25),(.65,1.,1.35,1.75)):
        configs.append(arm(f'joint_t{t}_s{s}_xy{xy}',t,s,xy,xy))
    return configs


def swept():
    return {k:sorted({float(a['env'][k]) for a in screen()}) for k in KEYS}


def extension(plan):
    configs = [arm('baseline')]
    spec = plan['extension']
    for t,s in itertools.product(spec['ridge_tau'],spec['ridge_sigma']):
        # Only a search coordinate. Cadence/regularizer clamps and nonlinear
        # state coupling prevent any claim of exact filter invariance.
        xy = round(1.15*t**(-41./14.),9)
        configs.append(arm(f'ext_ridge_t{t}_s{s}',t,s,xy,xy))
    for s,xy in itertools.product(spec['face_sigma'],spec['upper_xy']):
        configs.append(arm(f'ext_upper_s{s}_xy{xy}',.9,s,xy,xy))
    for s,xy in itertools.product(spec['face_sigma'],spec['lower_xy']):
        configs.append(arm(f'ext_lower_s{s}_xy{xy}',1.15,s,xy,xy))
    for s,xy in itertools.product(spec['coupled_sigma'],spec['coupled_xy']):
        configs.append(arm(f'ext_coupled_s{s}_xy{xy}',1.05,s,xy,xy))
    return configs


def split(plan):
    """Independent horizontal split, searched where OU-III's own split lives.

    TFG applies this knob exactly as OU-III does, set_RS_noise(rs*k_Sx, rs*k_Sy,
    rs), and still carries k_Sx = k_Sy = 1.15. Deployed OU-III carries 0.72 and
    0.50 on these same eight records, for a reason that belongs to the dataset
    rather than to that filter: yaw is pinned at 0 in every record, so world x
    and y are the vessel's surge and sway, and the dataset's keel damps sway but
    not surge. R_S is applied in world NED in both filters, so the reason
    transfers unchanged.

    The screen moves both horizontal axes together and so averages two axes that
    may disagree, and the refinement box reaches only 0.85. The pre-RAO note
    rejecting 0.72 for TFG tested x = y = 0.72; an isotropic rejection is not
    evidence about a split, which is the distinction OU-III had to draw when its
    own isotropic 0.5 sweep was rejected and the split was not. This grid
    therefore reaches well below the pre-RAO point on both axes independently.
    """
    spec = plan['split']
    configs = [arm('baseline')]
    for x,y in itertools.product(spec['TFG_R_S_X_FACTOR'],spec['TFG_R_S_Y_FACTOR']):
        if coeffs(arm('probe',spec['tau'],spec['sigma'],x,y)['env'])==LEGACY:
            continue
        configs.append(arm(f'split_x{x}_y{y}',spec['tau'],spec['sigma'],x,y))
    return configs


def refinement(plan):
    configs = [arm('baseline')]
    for t,s,x,y in itertools.product(*(plan['refinement'][k] for k in KEYS)):
        configs.append(arm(f'fine_t{t}_s{s}_x{x}_y{y}',t,s,x,y))
    return configs


def configs_for(phase,plan):
    configs = {'screen':screen,'extend':lambda:extension(plan),'split':lambda:split(plan),
               'refine':lambda:refinement(plan)}[phase]()
    values = [coeffs(c['env']) for c in configs]
    if len(set(values))!=len(values):
        raise ValueError('duplicate coefficient tuples in phase plan')
    return configs


def load(out,tags):
    rows = []
    seen = set()
    for tag in tags:
        path = out/tag/'runs.json'
        if not path.exists():
            continue
        for row in json.loads(path.read_text()):
            key = (row['config'],row['input'],str(row['seed']))
            if key not in seen:
                seen.add(key)
                rows.append(row)
    if not rows:
        raise ValueError('no completed requested phase')
    return rows


def main():
    from tfg_rao_refit_analysis import analyse
    from tfg_rao_refit_replay import execute,write_json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'reports/results/tfg_rao_core_refit')
    parser.add_argument('--jobs',type=int,default=4)
    parser.add_argument('--phase',choices=('screen','extend','split','refine','freeze','holdout'),default='screen')
    parser.add_argument('--plan',type=Path,default=ROOT/'tools/tfg_rao_refit_plan.json')
    parser.add_argument('--selection-runs',type=Path)
    parser.add_argument('--candidate-config')
    parser.add_argument('--cache-dir',type=Path)
    parser.add_argument('--shard-index',type=int,default=0)
    parser.add_argument('--shard-count',type=int,default=1)
    args = parser.parse_args()
    if not 0<=args.shard_index<args.shard_count:
        parser.error('invalid shard index/count')
    out = args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    cache = args.cache_dir or out/'cache'
    if args.phase in ('freeze','holdout'):
        if args.shard_count!=1:
            parser.error('freeze and holdout cannot be sharded by this entry point')
        frozen = out/'candidate.json'
        if args.phase=='freeze':
            if frozen.exists() or (out/'holdout-opened.json').exists():
                raise ValueError('candidate already frozen or holdout already opened')
            if args.selection_runs is None:
                parser.error('--selection-runs is required for a freeze')
            selection_bytes = args.selection_runs.read_bytes()
            selection = json.loads(selection_bytes)
            if any(str(r['seed']) in SEEDS['holdout'] for r in selection):
                raise ValueError('holdout data cannot be used for selection')
            analyse(selection)  # Require the complete eight-record product.
            ranking = rank(selection)
            chosen = next((r for r in ranking if r['eligible'] and
                           (args.candidate_config is None or r['config']==args.candidate_config)),None)
            if chosen is None:
                raise ValueError('no complete eligible requested candidate')
            coefficients = chosen['coefficients']
            payload = dict(config=chosen['config'],coefficients=coefficients,
                           coefficient_sha256=hashlib.sha256(json.dumps(coefficients,sort_keys=True).encode()).hexdigest(),
                           selection_runs_sha256=hashlib.sha256(selection_bytes).hexdigest(),
                           source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                           selected_on='training/refinement only',refinement_score=chosen['score'],
                           holdout_seeds_used_for_selection=False)
            with frozen.open('x') as handle:
                json.dump(payload,handle,indent=2)
                handle.write('\n')
            print(json.dumps(payload,indent=2))
            return
        candidate = json.loads(frozen.read_text())
        coefficients = candidate['coefficients']
        digest = hashlib.sha256(json.dumps(coefficients,sort_keys=True).encode()).hexdigest()
        if digest!=candidate['coefficient_sha256']:
            raise ValueError('frozen candidate checksum mismatch')
        opened = out/'holdout-opened.json'
        if opened.exists():
            if json.loads(opened.read_text())['coefficient_sha256']!=digest:
                raise ValueError('cannot change candidate on an opened holdout')
        else:
            with opened.open('x') as handle:
                json.dump({'coefficient_sha256':digest,'seeds':SEEDS['holdout']},handle)
        configs = [arm('baseline'),arm('candidate',*(coefficients[k] for k in KEYS))]
        rows = execute(configs,'TFG',PROTOCOL['record_indices'],SEEDS['holdout'],out/'holdout',args.jobs,cache)
        write_json(out/'holdout-analysis.json',analyse(rows))
        write_json(out/'holdout-ranking.json',rank(rows))
        execute([{'name':'ou3_shipping','env':{}}],'OU-III',PROTOCOL['record_indices'],
                SEEDS['holdout'],out/'holdout-ou3',args.jobs,cache)
        return
    plan = json.loads(args.plan.read_text()) if args.phase!='screen' else {}
    configs = configs_for(args.phase,plan)[args.shard_index::args.shard_count]
    if not configs:
        parser.error('empty shard')
    seeds = SEEDS[args.phase]
    if set(seeds)&set(SEEDS['holdout']):
        raise ValueError('selection and holdout seeds overlap')
    rows = execute(configs,'TFG',PROTOCOL['record_indices'],seeds,out/args.phase,args.jobs,cache)
    if args.shard_count==1:
        write_json(out/'analysis.json',analyse(rows))
        write_json(out/(args.phase+'-ranking.json'),rank(rows))
    write_json(out/'phase-status.json',dict(phase=args.phase,shard_index=args.shard_index,
               shard_count=args.shard_count,configurations=len(configs),records=8,seeds=seeds,
               holdout_evaluated=False,source_broad_run=plan.get('source_broad_run')))


if __name__=='__main__':
    main()
