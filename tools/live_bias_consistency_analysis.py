#!/usr/bin/env python3
"""Seed-level comparisons and bounded scientific figures; never scores plots."""
from pathlib import Path
import json,re,sys
import numpy as np
import pandas as pd
from scipy.stats import t as student
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/results/live_bias_consistency'
sys.path.insert(0,str(ROOT/'plots'))
import figure_decimation
from vector_export import save_figure
figure_decimation.install()
COLORS={'endpoint':'#596579','core':'#0a8292','frontend':'#c17b2d','both':'#6b529a','native':'#596579','schmidt':'#0a8292','schmidt_prior':'#6b529a'}

def ci(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v);m=float(v.mean()) if n else np.nan
 half=float(student.ppf(.975,n-1)*v.std(ddof=1)/np.sqrt(n)) if n>1 else 0.
 return dict(n=n,mean=m,low=m-half,high=m+half)
def label(v):
 x=ci(v);return f"{x['mean']:.2f} [{x['low']:.2f}, {x['high']:.2f}]"
def save(fig,name):
 fig.tight_layout();rows=[]
 def capture(_):
  rows.clear()
  for ai,ax in enumerate(fig.axes):
   for li,line in enumerate(ax.lines):
    x,y=line.get_data();rows.extend(dict(axis=ai,line=li,label=line.get_label(),x=a,y=b) for a,b in zip(x,y))
 cid=fig.canvas.mpl_connect('draw_event',capture)
 svg=OUT/'figures'/f'{name}.svg'
 save_figure(fig,svg,metadata={'Date':None})
 svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
 fig.canvas.mpl_disconnect(cid)
 pd.DataFrame(rows).to_csv(OUT/'figures'/f'{name}_vertices.csv.gz',index=False,float_format='%.9g',compression={'method':'gzip','mtime':0})
 save_figure(fig,Path('/tmp')/f'ocean_bias_{name}.png',dpi=140)
 plt.close(fig)
def read(stage):return pd.read_csv(OUT/f'{stage}_runs.csv')
def point(ax,x,values,color,label_=None):
 c=ci(values);ax.errorbar([x],[c['mean']],yerr=[[c['mean']-c['low']],[c['high']-c['mean']]],fmt='o',capsize=3,color=color,label=label_)

def clean_json(value):
 if isinstance(value,dict):return {key:clean_json(v) for key,v in value.items()}
 if isinstance(value,(list,tuple)):return [clean_json(v) for v in value]
 if isinstance(value,np.generic):value=value.item()
 return None if isinstance(value,float) and not np.isfinite(value) else value

def main():
 (OUT/'figures').mkdir(exist_ok=True)
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'live-bias-consistency'})
 stages={p.stem[:-5]:pd.read_csv(p) for p in OUT.glob('*_runs.csv')}
 run_keys=['seed','filter','motion','arm','timing','freeze','injection','duration','noise','lag']
 for name,frame in stages.items():
  manifest=json.loads((OUT/f'{name}_manifest.json').read_text())
  # The pre-latency stages predate that column and all used zero delay.
  keys=frame.assign(lag=frame['lag'] if 'lag' in frame else 0)
  assert len(keys[run_keys].drop_duplicates())==manifest['runs'],f'incomplete stage: {name}'
 assert all(k in stages for k in ('baseline','screen','long','recovery','noise'))
 b,s,d=stages['baseline'],stages['screen'],stages['long']
 late=d[(d.start==5400)|((d.motion=='rocking')&(d.start==2400))]
 blate=b[(b.start==5400)|((b.motion=='rocking')&(b.start==2400))]
 summary={};sections=[]
 sections+=['Six independent seeds; means [Student 95% intervals]. Bias units are mm/s².','', '| Condition | QMEKF | TFG |','| --- | ---: | ---: |']
 cases=[('Wave, latent acceleration, unchanged main','baseline','wave','exact_fixed','endpoint'),('Wave, latent acceleration, centered prediction','long','wave','exact_fixed','both'),('Wave, acceleration oracle, unchanged main','baseline','wave','oracle','endpoint'),('Wave, acceleration oracle, centered prediction','long','wave','oracle','both'),('Independent rocking, no oracle, centered prediction','long','rocking','exact_fixed','both')]
 values={}
 for title,stage,motion,arm,timing in cases:
  rows=[]
  for f in ('qmekf','tfg'):
   q=(blate if stage=='baseline' else late);v=q[(q['filter']==f)&(q.motion==motion)&(q.arm==arm)&(q.timing==timing)].ba_h_rms.to_numpy()*1000
   assert len(v)==6,(title,f,len(v));values[title,f]=v;rows.append(label(v));summary[title+' '+f]=ci(v)
  sections.append('| '+title+' | '+' | '.join(rows)+' |')
 sections+=['','Wave windows are 5400–7200 s; rocking is 2400–3600 s. Exact magnetic reference and fixed correct calibration are diagnostic controls. The oracle is artificial.','']
 effects=[]
 for f in ('qmekf','tfg'):
  q=d[(d['filter']==f)&(d.start==5400)&(d.arm=='exact_fixed')].pivot(index='seed',columns='timing',values='ba_h_rms')
  v=(q.both-q.endpoint)*1000;effects.append(f"{f.upper()}: {label(v)}")
  summary['paired timing effect '+f]=ci(v)
 sections+=['The paired centered-minus-endpoint effect with acceleration latent is '+ '; '.join(effects)+'. Both intervals include zero.','', '![Late bias errors](../reports/results/live_bias_consistency/figures/late_bias.svg)','']
 fig,axs=plt.subplots(1,2,figsize=(9,3.5),sharey=True)
 short=['Wave\nlatent','Wave\ncentered','Oracle\nendpoint','Oracle\ncentered','Rocking\ncentered']
 for ax,f in zip(axs,('qmekf','tfg')):
  for i,(title,*_) in enumerate(cases):point(ax,i,values[title,f],['#596579','#0a8292','#bd7b2a','#6b529a','#2a8b52'][i])
  ax.set_xticks(range(5),short);ax.set_title(f.upper());ax.grid(axis='y',alpha=.2)
 axs[0].set_ylabel('Horizontal bias RMS [mm/s²]');save(fig,'late_bias')
 # Localize core and frontend with common histories/tuning.
 screen=s[(s.start==900)&(s.motion=='wave')&(s.freeze=='none')]
 fig,axs=plt.subplots(1,2,figsize=(9,3.2))
 for ax,arm in zip(axs,('exact_fixed','oracle')):
  for fi,f in enumerate(('qmekf','tfg')):
   for i,timing in enumerate(('endpoint','core','frontend','both')):
    v=screen[(screen['filter']==f)&(screen.arm==arm)&(screen.timing==timing)].ba_h_rms*1000
    point(ax,i+(fi-.5)*.12,v,('#596579','#0a8292')[fi],f.upper() if i==0 else None)
  ax.set_xticks(range(4),['Endpoint','Core only','Frontend only','Both']);ax.set_title('Latent acceleration' if arm=='exact_fixed' else 'Acceleration oracle');ax.grid(axis='y',alpha=.2);ax.legend()
 axs[0].set_ylabel('900–1200 s horizontal bias RMS [mm/s²]');save(fig,'timing_localization')
 sections+=['![Core versus frontend localization](../reports/results/live_bias_consistency/figures/timing_localization.svg)','The core-only correction reproduces the oracle benefit; the frontend-only correction does not.','']
 # Independent freezes, paired by seed.
 fig,axs=plt.subplots(1,2,figsize=(8,3.2),sharey=True);freeze_stats=[]
 for ax,f in zip(axs,('qmekf','tfg')):
  q=s[(s.start==900)&(s.motion=='wave')&(s.arm=='exact_fixed')&(s.timing=='both')&(s['filter']==f)].pivot(index='seed',columns='freeze',values='ba_h_rms')
  for i,freeze in enumerate(('sea','coupling','both')):
   v=(q[freeze]-q.none)*1000;point(ax,i,v,'#0a8292');freeze_stats.append(dict(filter=f,freeze=freeze,**ci(v)))
  ax.axhline(0,color='#777',ls='--',lw=.8);ax.set_xticks(range(3),['Sea-state hold','Coupling hold','Both holds']);ax.set_title(f.upper());ax.grid(axis='y',alpha=.2)
 axs[0].set_ylabel('Paired change in horizontal bias RMS [mm/s²]');save(fig,'adaptation_freezes');summary['freeze_effects']=freeze_stats
 sections+=['![Independent adaptation freezes](../reports/results/live_bias_consistency/figures/adaptation_freezes.svg)','']
 # Seed 1001 trace: descriptive, never used to score RMS.
 fig,axs=plt.subplots(2,2,figsize=(9,5.5),sharex=True)
 for col,f in enumerate(('qmekf','tfg')):
  for arm,timing,style in [('exact_fixed','endpoint','--'),('exact_fixed','both','-'),('oracle','endpoint','--'),('oracle','both','-')]:
   p=OUT/'runs/long'/f'1001_{f}_wave_{arm}_{timing}_none_none_7200_trace.csv';q=pd.read_csv(p)
   e=np.hypot(q.ba_x-q.true_x,q.ba_y-q.true_y)*1000;color='#596579' if arm=='exact_fixed' else '#0a8292'
   title=('Latent' if arm=='exact_fixed' else 'Oracle')+' / '+timing
   axs[0,col].plot(q.t/60,e,style,color=color,lw=.9,label=title);axs[1,col].plot(q.t/60,q.tilt,style,color=color,lw=.9)
  axs[0,col].set_title(f.upper()+' — seed 1001');axs[0,col].legend(fontsize=7);axs[1,col].set_xlabel('Time [min]')
  for ax in axs[:,col]:ax.axvline(10,color='#aaa',lw=.7);ax.grid(alpha=.2)
 axs[0,0].set_ylabel('Horizontal bias error [mm/s²]');axs[1,0].set_ylabel('Tilt error [deg]');save(fig,'long_traces')
 sections+=['![Representative two-hour trace](../reports/results/live_bias_consistency/figures/long_traces.svg)','']
 # Coverage/NEES and properly tick-scored innovations.
 consistency=[]
 for (f,arm,timing),q in late[late.motion=='wave'].groupby(['filter','arm','timing']):
  consistency.append(dict(filter=f,arm=arm,timing=timing,coverage=ci(q.coverage),h_nees=ci(q.h_nees),acc_nis=ci(q.acc_nis),mag_nis=ci(q.mag_nis),enabled=ci(q.enabled_fraction),applied=ci(q.mag_applied_fraction)))
 summary['consistency']=consistency
 sections+=['Covariance diagnostics remain conditional on the model. Late oracle coverage changes from 0% to 100% in these serial samples when timing is corrected; the latent-wave errors remain inconsistently covered. A 100% serial coverage fraction is not a calibration guarantee. Actual bias gates are open in scored live windows. Magnetic NIS is accumulated only at packets, and applied-correction counters are recorded separately.','']
 # Window information and relinearization, seeds remain independent units.
 windows=json.loads((OUT/'window_results.json').read_text());assert len(windows)==72
 wd=pd.DataFrame([dict(seed=x['seed'],filter=x['filter'],motion=x['motion'],case=x['case'],initial=x['initial']['ba_h_mm_s2'],one_step=x['one_step']['ba_h_mm_s2'],relinearized=x['relinearized']['ba_h_mm_s2'],gradient=x['optimality'],status=x['status'],weak_information=x['tilt_bias_motion_information'][1]['eigenvalues'][0]) for x in windows])
 wd.to_csv(OUT/'window_statistics.csv',index=False,float_format='%.10g')
 sections+=['| 60-s window, QMEKF initialization | Live initialization | First linearized solve | Relinearized |','| --- | ---: | ---: | ---: |']
 for motion in ('wave','rocking'):
  for case in ('native_prior','motion4','bias3'):
   q=wd[(wd['filter']=='qmekf')&(wd.motion==motion)&(wd['case']==case)];sections.append(f'| {motion}; {case} | {label(q.initial)} | {label(q.one_step)} | {label(q.relinearized)} |')
 sections+=['','TFG initializations give nearly identical window solutions. The window uses the known synthetic white-noise parameters and its declared population bias prior. A change from live initialization includes that changed window objective; only first-solve versus converged results isolate further relinearization. All 72 fits converge; PSD floors are unused. Weak information and its regularizer sensitivity are retained in JSON. The five-dimensional tilt/bias analysis eliminates world yaw as well as all other poses, gyro bias and motion. Prior sensitivity and degraded rocking recovery do not qualify window feedback.','']
 fig,axs=plt.subplots(1,2,figsize=(8,3.2))
 for ax,motion in zip(axs,('wave','rocking')):
  for i,case in enumerate(('native_prior','motion4','bias3')):
   q=wd[(wd['filter']=='qmekf')&(wd.motion==motion)&(wd['case']==case)];point(ax,i-.12,q.initial,'#596579','Live initialization' if i==0 else None);point(ax,i,q.one_step,'#bd7b2a','First solve' if i==0 else None);point(ax,i+.12,q.relinearized,'#0a8292','Relinearized' if i==0 else None)
  ax.set_title(motion.capitalize());ax.set_xticks(range(3),['Declared priors','4× motion cov.','3× bias SD']);ax.grid(axis='y',alpha=.2);ax.legend(fontsize=7)
 axs[0].set_ylabel('Window horizontal bias error [mm/s²]');save(fig,'window_relinearization')
 sections+=['![Window prior sensitivity](../reports/results/live_bias_consistency/figures/window_relinearization.svg)','']
 fig,axs=plt.subplots(1,2,figsize=(8,3.2),sharey=True)
 for ax,motion in zip(axs,('wave','rocking')):
  selected=[w for w in windows if w['filter']=='qmekf' and w['motion']==motion and w['case']=='native_prior']
  for key,title,color in [('tilt_bias_measurement_information','Sensors only','#596579'),('tilt_bias_motion_information','+ motion prior','#0a8292'),('tilt_bias_information','+ bias prior','#6b529a')]:
   eig=np.median([w[key][1]['eigenvalues'] for w in selected],axis=0)
   ax.plot(range(1,6),np.maximum(eig,1e-10),'o-',color=color,label=title)
  ax.set_title(motion.capitalize()+' — six-seed median');ax.set_yscale('log');ax.set_xticks(range(1,6));ax.set_xlabel('Ordered eigenvalue');ax.grid(alpha=.2);ax.legend(fontsize=7)
 axs[0].set_ylabel('Scaled tilt/bias information (yaw eliminated)');save(fig,'window_information')
 sections+=['![Window information and prior contribution](../reports/results/live_bias_consistency/figures/window_information.svg)','Eigenvalues at or below 1e-10 are plotted at that numerical floor. The sensor-only model in this plot deliberately frees the motion history.','']
 # Optional later stages make the script useful while a long experiment runs.
 for stage,start in [('noise',900),('calibration',2400),('combined',5400),('latency',900),('startup',900),('recovery',2400),('late_recovery',5400)]:
  if stage not in stages:continue
  q=stages[stage];q=q[q.start==start]
  groups=['filter','motion','arm','timing']+[c for c in ('noise','lag','injection') if c in q]
  stats=[]
  for key,g in q.groupby(groups):
   row=dict(zip(groups,key));row.update(bias=ci(g.ba_h_rms*1000),tilt=ci(g.tilt_rms),h_nees=ci(g.h_nees),coverage=ci(g.coverage),mag_nis=ci(g.mag_nis));stats.append(row)
  summary[stage]=stats
 if 'calibration' in stages:
  q=stages['calibration'];q=q[(q.start==2400)&(q.arm=='exact_fixed')];fig,axs=plt.subplots(1,2,figsize=(9,3.2))
  kinds=['prior','late_prior','all_prior','schmidt','schmidt_prior']
  for ax,f in zip(axs,('qmekf','tfg')):
   for i,kind in enumerate(kinds):point(ax,i,q[(q['filter']==f)&(q.noise==kind)].ba_h_rms*1000,'#0a8292')
   ax.set_title(f.upper());ax.set_xticks(range(5),['Initial prior','Late prior','White + prior','Schmidt','Schmidt + prior'],rotation=18);ax.grid(axis='y',alpha=.2)
  axs[0].set_ylabel('2400–3600 s bias RMS [mm/s²]');save(fig,'calibration_controls')
  sections+=['![Calibration uncertainty controls](../reports/results/live_bias_consistency/figures/calibration_controls.svg)','']
 if 'combined' in stages:
  q=stages['combined'];q=q[q.start==5400];sections+=['| Timing from startup; two-hour tail | QMEKF | TFG |','| --- | ---: | ---: |']
  for arm in ('exact_fixed','learned_fixed'):
   for noise in ('native','schmidt','schmidt_prior'):
    vals=[]
    for f in ('qmekf','tfg'):vals.append(label(q[(q['filter']==f)&(q.arm==arm)&(q.noise==noise)].ba_h_rms*1000))
    sections.append('| '+arm+'; '+noise+' | '+' | '.join(vals)+' |')
  sections+=['','These rows distinguish actual bias recovery from corrected covariance handling; the learned-reference rows use no truth for estimation.','']
 if 'production' in stages:
  q=stages['production'];tail=q[(q.start==5400)|((q.motion=='rocking')&(q.start==2400))]
  sections+=['| Actual default configuration | QMEKF | TFG |','| --- | ---: | ---: |']
  production_stats=[]
  for motion,arm in [('wave','exact_fixed'),('wave','learned_fixed'),('wave','oracle'),('rocking','exact_fixed')]:
   vals=[]
   for f in ('qmekf','tfg'):
    g=tail[(tail['filter']==f)&(tail.motion==motion)&(tail.arm==arm)]
    assert len(g)==6
    vals.append(label(g.ba_h_rms*1000))
    production_stats.append(dict(filter=f,motion=motion,arm=arm,bias=ci(g.ba_h_rms*1000),tilt=ci(g.tilt_rms),h_nees=ci(g.h_nees),coverage=ci(g.coverage),enabled=ci(g.enabled_fraction),mag_applied=ci(g.mag_applied_fraction)))
   sections.append('| '+motion+'; '+arm+' | '+' | '.join(vals)+' |')
  summary['production_defaults']=production_stats
  sections+=['','These runs leave sampling and held-bias Config defaults untouched under the common benchmark tuning: centered live prediction, the existing endpoint frontend, and legacy held-bias handling. The focused 180-s samplewise Config/setter control checks full state and covariance through startup and activation. Only the oracle arm supplies acceleration truth. The corrected implementation still has a large unknown-acceleration floor.','']
 if (OUT/'main_rollout_scores.csv').exists():
  q=pd.read_csv(OUT/'main_rollout_scores.csv',dtype={'sensor_hash':str})
  manifest=json.loads((OUT/'main_rollout_manifest.json').read_text())
  assert len(q)==manifest['filter_runs']==36 and (q.groupby(['family','seed']).sensor_hash.nunique()==1).all()
  assert manifest['duration']==7200 and (q.duration==7200).all() and (q.start==5400).all(),'rollout requires the full 5400–7200 s scoring window'
  sections+=['| Other live filters: acquired reference, unknown acceleration | Legacy RMS | Default RMS | Paired default minus legacy |','| --- | ---: | ---: | ---: |']
  rollout=[];fig,axs=plt.subplots(1,2,figsize=(9,3.4))
  for i,f in enumerate(('OU2','OU3','TFG')):
   g=q[q.family==f];p=g.pivot(index='seed',columns='profile',values='ba_h_rms')*1000
   sections.append('| '+f+' | '+label(p.legacy)+' | '+label(p.default)+' | '+label(p.default-p.legacy)+' |')
   row=dict(family=f,bias_legacy=ci(p.legacy),bias_default=ci(p.default),paired_bias_change=ci(p.default-p.legacy))
   for metric in ('tilt_rms','bg_rms','h_nees','coverage','acc_nis','enabled_fraction','applied_mag'):
    row[metric]={profile:ci(g[g.profile==profile][metric]) for profile in ('legacy','default')}
   rollout.append(row)
   for j,profile in enumerate(('legacy','default')):
    point(axs[0],i+(j-.5)*.15,p[profile],('#596579','#0a8292')[j],profile.capitalize() if i==0 else None)
    point(axs[1],i+(j-.5)*.15,g[g.profile==profile].tilt_rms,('#596579','#0a8292')[j],profile.capitalize() if i==0 else None)
  for ax in axs:ax.set_xticks(range(3),['OU-II','OU-III','TFG']);ax.grid(axis='y',alpha=.2);ax.legend()
  axs[0].set_ylabel('Horizontal bias RMS [mm/s²]');axs[1].set_ylabel('Tilt RMS [deg]');save(fig,'main_filter_rollout')
  summary['main_filter_rollout']=rollout
  sections+=['','Six paired seeds, 7200 s each, scoring 5400–7200 s. Each legacy/default pair consumes identical sensors from one continuous physical history. Magnetic reference is acquired ordinarily; fixed hard iron is correct. These are within-family timing comparisons with existing family priors, not cross-family rankings or exact-reference CARMA controls. The timing correction slightly worsens this particular floor in all three families. Covariance remains overconfident (zero late 95% horizontal coverage here), and actual bias learning is enabled throughout. Correct integration does not create missing independent motion/reference information.','', '![Other live filter rollout](../reports/results/live_bias_consistency/figures/main_filter_rollout.svg)','']
 if (OUT/'csv_rollout_all.json').exists():
  records=json.loads((OUT/'csv_rollout_all.json').read_text());csvstats=[]
  for r in records:
   source=ROOT/'tests'/('kalman_carma_shared_modes/kalman_carma_shared_modes-sim.cpp' if r['filter']=='qmekf' else 'kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim.cpp')
   limits={k:float(v) for k,v in re.findall(r'\.(\w+)\s*=\s*([0-9.]+)f',source.read_text())}
   line=next(x for x in r['metrics'] if x.startswith('VALIDATION_METRICS '));v=dict(x.split('=',1) for x in line.split()[1:])
   acc=next(x for x in r['metrics'] if x.startswith('Bias error RMS (% of max TRUE bias) (acc)'))
   gyro=next(x for x in r['metrics'] if x.startswith('Bias error RMS (% of max TRUE bias) (gyro)'))
   kind='jonswap' if 'jonswap' in r['input'] else 'pmstokes'
   metrics={'err_limit_percent_z_'+kind:float(v['disp_z_pct_hs']),'err_limit_percent_3d_'+kind:float(v['disp_3d_pct_refmax']),**{'err_limit_'+k+'_deg':float(v[k+'_rms_deg']) for k in ('roll','pitch','yaw')},'acc_z_bias_percent':float(re.search(r'Z=([0-9.]+)%',acc)[1]),'bias_3d_percent':float(re.search(r'\|3D\|=([0-9.]+)%',acc)[1]),'gyro_bias_3d_percent':float(re.search(r'\|3D\|=([0-9.]+)%',gyro)[1])}
   failures=[k for k,x in metrics.items() if x>limits[k]]
   csvstats.append(dict(filter=r['filter'],profile=r['contract'],input=r['input'],exit_code=r['exit_code'],failed_limits=';'.join(failures),roll=float(v['roll_rms_deg']),pitch=float(v['pitch_rms_deg']),yaw=float(v['yaw_rms_deg']),bias_mm=float(v['accel_bias_3d_rms_mps2'])*1000,gyro_bias=float(v['gyro_bias_3d_rms_radps'])))
  cs=pd.DataFrame(csvstats);cs.to_csv(OUT/'csv_rollout_statistics.csv',index=False,float_format='%.10g')
  sections+=['| Pinned replay profile | QMEKF records exceeding old bars | TFG records exceeding old bars |','| --- | ---: | ---: |']
  summaries=[]
  for profile in ('interval','defaults','core','frontend','instantaneous','core_schmidt','both_schmidt'):
   if profile not in ('interval','defaults') and not (cs.profile==profile).any():continue
   counts=[]
   for f in ('qmekf','tfg'):
    g=cs[(cs['filter']==f)&(cs.profile==profile)];assert len(g)==8
    counts.append(str(int((g.failed_limits!='').sum()))+'/8')
    summaries.append(dict(filter=f,profile=profile,failed_records=g[g.failed_limits!=''].input.to_list(),max_roll=float(g.roll.max()),max_pitch=float(g.pitch.max()),max_yaw=float(g.yaw.max()),max_bias_mm=float(g.bias_mm.max())))
    if profile in ('interval','defaults'):assert (g.exit_code==0).all() and (g.failed_limits=='').all(),'gated replay regressed'
   sections.append('| '+profile+' | '+' | '.join(counts)+' |')
  summary['csv_rollout']=summaries
  sections+=['','Each profile uses the same eight pinned histories and default noise seeds. `interval` is the locked baseline; `defaults` retains the actual sampling and held-bias Config defaults; `core` selects that sampling path explicitly. `instantaneous` centers both paths, and the two Schmidt profiles add held-bias handling. Diagnostic profiles are compared offline to the original bars, without claiming their ungated exit status certifies those bars. Production and legacy profiles are gated natively. This evidence selects core-only prediction and retains the original frontend/holding defaults. It does not establish the exact mechanism of the frontend/reference regression.','']
 # Covariance diagnostics use serial samples and their actual marginal P.
 fig,axs=plt.subplots(2,2,figsize=(9,5),sharex=True)
 for col,f in enumerate(('qmekf','tfg')):
  for arm,timing,style in [('exact_fixed','endpoint','--'),('exact_fixed','both','-'),('oracle','endpoint','--'),('oracle','both','-')]:
   q=pd.read_csv(OUT/'runs/long'/f'1001_{f}_wave_{arm}_{timing}_none_none_7200_diag.csv')
   color='#596579' if arm=='exact_fixed' else '#0a8292';name=('Latent' if arm=='exact_fixed' else 'Oracle')+' / '+timing
   axs[0,col].plot(q.t/60,q.h_nees,style,color=color,label=name,lw=.9)
   axs[1,col].plot(q.t/60,q.nuisance_inflation,style,color=color,lw=.9)
  axs[0,col].axhline(5.99146,color='#bd7b2a',lw=.8,label='95% ellipse threshold');axs[0,col].set_yscale('log');axs[0,col].set_title(f.upper()+' — seed 1001');axs[0,col].legend(fontsize=6)
  axs[1,col].set_xlabel('Time [min]')
  for ax in axs[:,col]:ax.set_xlim(2.5,120);ax.grid(alpha=.2)
 axs[0,0].set_ylabel('Horizontal bias NEES');axs[1,0].set_ylabel('Marginal / conditional variance ratio');save(fig,'covariance_diagnostics')
 sections+=['![Covariance consistency and nuisance conditioning](../reports/results/live_bias_consistency/figures/covariance_diagnostics.svg)','The variance ratio is a posterior-covariance diagnostic, not an observability Gramian.','']
 if 'noise' in stages:
  fig,axs=plt.subplots(1,2,figsize=(9,3.2))
  q=stages['noise'];q=q[(q.start==900)&(q.motion=='wave')&(q.timing=='both')]
  for ax,metric in zip(axs,('acc_nis','mag_nis')):
   for fi,f in enumerate(('qmekf','tfg')):
    for i,kind in enumerate(('native','gyro','mag','all')):point(ax,i+(fi-.5)*.12,q[(q['filter']==f)&(q.noise==kind)][metric],('#596579','#0a8292')[fi],f.upper() if i==0 else None)
   ax.axhline(3,color='#bd7b2a',lw=.8,label='Dimension reference');ax.set_xticks(range(4),['Native weights','Gyro matched','Mag matched','All matched']);ax.set_title(metric.replace('_',' ').upper());ax.grid(axis='y',alpha=.2);ax.legend(fontsize=7)
  axs[0].set_ylabel('900–1200 s mean NIS');save(fig,'innovations')
  sections+=['![Innovations under independently changed white-noise models](../reports/results/live_bias_consistency/figures/innovations.svg)','Matching magnetic white-noise SD can bring magnetic NIS near its dimension reference without resolving bias error. The prior/reference and motion assumptions remain separate.','']
 for recovery_stage in ('late_recovery','production_recovery'):
  if recovery_stage not in stages:continue
  controls=stages[recovery_stage].query('injection=="none" and arm in ["exact_fixed","learned_fixed"]')
  reference_stage='combined' if recovery_stage=='late_recovery' else 'production'
  if reference_stage in stages:
   reference=stages['combined'].query('noise=="schmidt"') if recovery_stage=='late_recovery' else stages['production']
   common_controls=controls.merge(reference,on=['seed','filter','motion','arm','duration','start','end'],suffixes=('_recovery','_combined'))
   control_difference=float(np.max(abs(common_controls.ba_h_rms_recovery-common_controls.ba_h_rms_combined)))
   assert len(common_controls)==120 and control_difference==0,'startup configuration/setter controls differ'
   summary['startup_configuration_crosscheck' if recovery_stage=='late_recovery' else 'production_recovery_control']=dict(windows=len(common_controls),max_bias_rms_difference=control_difference)
  q=stages[recovery_stage];late=q[q.start==5400];sections+=['**Recovery profile: '+recovery_stage+'.**','', '| One-hour shock; 5400–7200 s bias RMS | QMEKF | TFG |','| --- | ---: | ---: |']
  for motion,arm,inj in [('wave','exact_fixed','none'),('wave','exact_fixed','bias'),('wave','exact_fixed','tilt'),('rocking','exact_fixed','bias'),('rocking','exact_fixed','tilt'),('wave','learned_fixed','none'),('wave','correct_reference','none')]:
   vals=[]
   for f in ('qmekf','tfg'):
    g=late[(late['filter']==f)&(late.motion==motion)&(late.arm==arm)&(late.injection==inj)];vals.append(label(g.ba_h_rms*1000))
   sections.append('| '+motion+'; '+(inj if inj!='none' else arm)+' | '+' | '.join(vals)+' |')
  profile='centered core/frontend and Schmidt holding' if recovery_stage=='late_recovery' else 'the actual core-only production defaults'
  sections+=['','All these runs use '+profile+' from startup with the common native prior. At 3600 s, bias/tilt shocks are additive estimated-state perturbations; the reference-repair arm instead changes only the magnetic reference. No live covariance reset accompanies an intervention.','']
  recovery=[]
  for key,g in late.groupby(['filter','motion','arm','injection']):
   if key[3]=='none' and key[2]!='correct_reference':continue
   crossing=pd.to_numeric(g.recovery_bias,errors='coerce')
   recovery.append(dict(zip(['filter','motion','arm','injection'],key))|dict(first_qualifying_count=int(crossing.notna().sum()),late_below_10mm=int((g.ba_h_rms<.01).sum()),first_qualifying_seconds=ci(crossing)))
  summary[recovery_stage+'_counts']=recovery
  fig,axs=plt.subplots(2,2,figsize=(9,5.3),sharex=True)
  for col,f in enumerate(('qmekf','tfg')):
   for row,motion in enumerate(('wave','rocking')):
    timing,noise=('both_start','_schmidt') if recovery_stage=='late_recovery' else ('defaults','')
    for inj,color in [('none','#596579'),('bias','#bd7b2a'),('tilt','#0a8292')]:
     path=OUT/'runs'/recovery_stage/f'1001_{f}_{motion}_exact_fixed_{timing}_none_{inj}_7200{noise}_trace.csv';trace=pd.read_csv(path)
     error=np.hypot(trace.ba_x-trace.true_x,trace.ba_y-trace.true_y)*1000
     axs[row,col].plot(trace.t/60,error,color=color,lw=.9,label={'none':'Control','bias':'Bias shock','tilt':'Tilt shock'}[inj])
    axs[row,col].axvline(60,color='#999',ls='--',lw=.8);axs[row,col].axhline(10,color='#6b529a',ls=':',lw=.8);axs[row,col].set_xlim(55,120);axs[row,col].set_yscale('log');axs[row,col].grid(alpha=.2);axs[row,col].set_title(f.upper()+' / '+motion+' — seed 1001');axs[row,col].legend(fontsize=7)
  axs[0,0].set_ylabel('Horizontal bias error [mm/s²]');axs[1,0].set_ylabel('Horizontal bias error [mm/s²]')
  for ax in axs[1]:ax.set_xlabel('Time [min]')
  save(fig,recovery_stage)
  sections+=['![One-hour shock recovery](../reports/results/live_bias_consistency/figures/'+recovery_stage+'.svg)','Recovery counts distinguish a first qualifying trailing 30-s window, sustained for another 30 s, from the independent late-window threshold. They are retained by seed in the summary.','']
 # Startup trace includes actual magnetic packet events near handoff.
 if 'startup' in stages:
  path=OUT/'runs/production'/'1001_qmekf_wave_learned_fixed_defaults_none_none_7200_trace.csv'
  if not path.exists():path=OUT/'runs/startup'/'1001_qmekf_wave_learned_fixed_both_none_none_1200_trace.csv'
  q=pd.read_csv(path);q=q[q.t<=180]
  fig,axs=plt.subplots(2,2,figsize=(9,5),sharex=True)
  axs[0,0].plot(q.t,q.tilt,color='#0a8292');axs[0,0].set_ylabel('Tilt error [deg]')
  axs[1,0].plot(q.t,np.hypot(q.ba_x-q.true_x,q.ba_y-q.true_y)*1000,color='#0a8292');axs[1,0].set_ylabel('Bias error [mm/s²]')
  for key,title,color in [('live','Live','#596579'),('refined','Refined reference','#bd7b2a'),('enabled','Bias gate','#0a8292')]:axs[0,1].step(q.t,q[key],where='post',label=title,color=color)
  axs[0,1].legend(fontsize=7);axs[0,1].set_ylabel('Actual flag')
  # The magnitude is retained between packets; accepted-event markers are
  # shown only where the native magnetic kernel actually ran.
  axs[1,1].plot(q.t,q.mag_correction_rad*180/np.pi,color='#596579',lw=.8,label='Most recent correction')
  applied=q[q.mag_applied.astype(bool)];axs[1,1].plot(applied.t,applied.mag_correction_rad*180/np.pi,'.',ms=2,color='#0a8292',label='Applied magnetic tick')
  axs[1,1].set_ylabel('Magnetic attitude correction [deg]');axs[1,1].legend(fontsize=7)
  for ax in axs.flat:ax.grid(alpha=.2)
  for ax in axs[1]:ax.set_xlabel('Time [s]')
  save(fig,'startup_gates')
  sections+=['![Startup and actual learning gates](../reports/results/live_bias_consistency/figures/startup_gates.svg)','']
 # Validate unchanged default on the long common prefixes/tails.
 assert 'compatibility' in stages,'compatibility_runs.csv is required for baseline preservation'
 comparison=stages['compatibility']
 common=comparison.merge(b,on=['seed','filter','motion','arm','start','end'],suffixes=('_new','_base'))
 summary['default_preservation']=dict(rows=len(common),max_bias_rms_difference=float(np.max(abs(common.ba_h_rms_new-common.ba_h_rms_base))))
 assert len(common)==168 and summary['default_preservation']['max_bias_rms_difference']==0,'default regression changed'
 sections+=['The untouched-main binary and the modified binary in compatibility mode match horizontal-bias RMS exactly on all '+str(len(common))+' overlapping scored windows. The baseline was separately compiled against untouched production headers.','']
 (OUT/'summary.json').write_text(json.dumps(clean_json(summary),indent=2,allow_nan=False)+'\n')
 doc=ROOT/'docs/live-bias-consistency.md';text=doc.read_text();a=text.index('<!-- RESULTS_START -->')+len('<!-- RESULTS_START -->');z=text.index('<!-- RESULTS_END -->');doc.write_text(text[:a]+'\n\n'+'\n'.join(sections)+'\n'+text[z:])
 print('Wrote report, seed-level statistics and bounded figures.')
if __name__=='__main__':main()
