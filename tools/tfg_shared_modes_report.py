#!/usr/bin/env python3
"""Bounded SVG/PGF figures, scientific tables and article macros for the TFG study."""
from __future__ import annotations
import argparse,json,re,sys,os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'plots'))
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.family':'serif','font.size':8,'axes.labelsize':8,
 'legend.fontsize':6.5,'xtick.labelsize':7,'ytick.labelsize':7,'axes.grid':True,
 'grid.alpha':.23,'axes.unicode_minus':False,'pgf.texsystem':'pdflatex','pgf.rcfonts':False})
import matplotlib.pyplot as plt
from vector_export import save_figure
from svg_determinism import configure_svg
configure_svg(matplotlib)
LABEL={'A41@C':'QMEKF independent','T41@C':'TFG independent','CSM':'QMEKF shared','TCSM':'TFG shared'}
COLOR=dict(zip(LABEL,['#777777','#009e73','#0072b2','#d55e00']))
CASE={'jonswap':'JONSWAP','pm_long':'PM, long period','regular':'Regular waves','short_low':'Short, low waves',
 'high_long':'High, long waves','broad_direction':'Broad spreading','long_crested':'Long crested',
 'heading90':'Heading +90 deg','head_seas':'Head seas','swell':'Directional swell','mixed':'Mixed seas',
 'weak_heave':'Reduced heave','calm':'Calm water','low_energy':'Low wave energy',
 'engine800':'Engine 800 rpm','engine1600':'Engine 1600 rpm','engine2400':'Engine 2400 rpm','engine3200':'Engine 3200 rpm',
 'colored_noise':'Colored sensor noise','roundtrip':'Sea-state round trip','direction_step':'Direction step',
 'changing_period':'Changing period','maneuver':'Vessel maneuvers','rapid_rotation':'Rapid rotations',
 'mag_outage':'Magnetic outage','bias_temperature':'Bias / temperature drift','jitter':'Interval jitter'}
OUT=ROOT/'reports/results/tfg_shared_modes';DOC=ROOT/'doc/kalman_tfg_shared_modes'
def tex(s):
 """Escape plain text used in generated LaTeX cells."""
 return str(s).replace('_',r'\_').replace('%',r'\%').replace('&',r'\&')
def fmt(x,n=2):
 """Format finite measurements, marking missing values explicitly."""
 return '--' if not np.isfinite(x) else f'{x:.{n}f}'
def interval(r,percent=True):
 """Render a paired percent change with its marginal confidence interval."""
 return f"{r['percent']:+.2f} [{r['pct_ci_low']:+.2f}, {r['pct_ci_high']:+.2f}]" if percent else ''
def save(fig,out,name):
 """Export the same bounded figure geometry to deterministic SVG and PGF."""
 fig.align_ylabels();fig.tight_layout(pad=.65)
 for ext in ('svg','pgf'):
  path=out/f'{name}.{ext}'
  save_figure(fig,path,bbox_inches='tight',**({'metadata':{'Date':None}} if ext=='svg' else {}))
  path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
 plt.close(fig)
def trace(results,name):
 """Load a decimated publication trace, never a scoring input."""
 return pd.read_csv(results/'traces'/f'{name}.csv.gz')
def paired_figures(results,out,pairs,interactions):
 """Plot paired position changes and geometry/prior interactions by condition."""
 cases=list(CASE);d=pairs[(pairs.comparison=='TFG geometry, shared')&pairs.metric.isin(['rms_z','rms_h'])]
 fig,axs=plt.subplots(1,2,figsize=(7.1,5.15),sharey=True)
 for ax,m,title in zip(axs,['rms_z','rms_h'],['Heave RMS','Horizontal RMS']):
  s=d[d.metric==m].set_index('case').reindex(cases);y=np.arange(len(cases));ax.axvline(0,color='black',lw=.7)
  ax.axvspan(-5,5,color='#cccccc',alpha=.18)
  ax.errorbar(s.percent,y,xerr=[s.percent-s.pct_ci_low,s.pct_ci_high-s.percent],fmt='o',ms=3,lw=.8,color=COLOR['TCSM'])
  ax.set_yticks(y,[CASE[x] for x in cases]);ax.set_title(title);ax.set_xlabel('TFG / QMEKF shared change [%]');ax.set_ylim(len(cases)-.5,-.5)
 save(fig,out,'paired_position')
 fig,axs=plt.subplots(1,2,figsize=(7.1,5.15),sharey=True)
 for ax,m,title in zip(axs,['rms_z','rms_h'],['Heave interaction','Horizontal interaction']):
  s=interactions[interactions.metric==m].set_index('case').reindex(cases);y=np.arange(len(cases));ax.axvline(0,color='black',lw=.7)
  ax.errorbar(1000*s.interaction,y,xerr=[1000*(s.interaction-s.ci_low),1000*(s.ci_high-s.interaction)],fmt='o',ms=3,lw=.8,color='#cc79a7')
  ax.set_yticks(y,[CASE[x] for x in cases]);ax.set_title(title);ax.set_xlabel('Difference in differences [mm]');ax.set_ylim(len(cases)-.5,-.5)
 save(fig,out,'factorial_interaction')
def tracking(results,out):
 """Plot final position, velocity, attitude and bias tracking samples."""
 d=trace(results,'jonswap_tail_position');fig,axs=plt.subplots(3,1,figsize=(3.45,4.3),sharex=True)
 for ax,a in zip(axs,'xyz'):
  ax.plot(d.t,d[a+'_true'],color='black',lw=.9,label='Truth')
  for f in LABEL:ax.plot(d.t,d[a+'_'+f],lw=.7,color=COLOR[f],ls='--' if f!='TCSM' else '-',label=LABEL[f])
  ax.set_ylabel(a.upper()+' [m]')
 axs[0].legend(ncol=2);axs[-1].set_xlabel('Time [s]');save(fig,out,'tail_xyz')
 d=trace(results,'jonswap_tail_velocity');fig,axs=plt.subplots(3,1,figsize=(3.45,4.0),sharex=True)
 for ax,a in zip(axs,'xyz'):
  ax.plot(d.t,d['v'+a+'_true'],color='black',lw=.9,label='Truth')
  for f in LABEL:ax.plot(d.t,d['v'+a+'_'+f],lw=.7,color=COLOR[f],label=LABEL[f])
  ax.set_ylabel(a.upper()+' velocity [m/s]')
 axs[0].legend(ncol=2);axs[-1].set_xlabel('Time [s]');save(fig,out,'tail_velocity')
 d=trace(results,'jonswap_tail_attitude');fig,axs=plt.subplots(3,1,figsize=(3.45,4.0),sharex=True)
 for ax,a in zip(axs,['roll','pitch','yaw']):
  ax.plot(d.t,d[a+'_true'],color='black',lw=.9,label='Truth')
  for f in LABEL:
   e=d[a+'_true']+d['d'+a+'_'+f];ax.plot(d.t,e,lw=.7,color=COLOR[f],label=LABEL[f])
  ax.set_ylabel(a.capitalize()+' [deg]')
 axs[0].legend(ncol=2);axs[-1].set_xlabel('Time [s]');save(fig,out,'attitude')
 for family,stem,scale,unit,name in [('gyro_bias','bg',180/np.pi*3600,'deg/h','Gyro'),('accel_bias','ba',1000/9.80665,'mg','Accel.')]:
  d=trace(results,'jonswap_'+family);fig,axs=plt.subplots(3,1,figsize=(3.45,4.0),sharex=True)
  for ax,a in zip(axs,'xyz'):
   ax.plot(d.t,scale*d[stem+a+'_true'],color='black',lw=.9,label='Truth')
   for f in LABEL:ax.plot(d.t,scale*d[stem+a+'_'+f],lw=.7,color=COLOR[f],label=LABEL[f])
   ax.set_ylabel(f'{name} {a.upper()} [{unit}]')
  axs[0].legend(ncol=2);axs[-1].set_xlabel('Time [s]');save(fig,out,family)
def degradation(results,out,df):
 """Plot engine and sea-transition degradation using the registered scores."""
 fig,axs=plt.subplots(2,1,figsize=(3.45,3.8),sharex=True)
 for f in LABEL:
  for ax,m in zip(axs,['rms_z','rms_h']):
   ys=[];err=[]
   for case in ['jonswap','engine800','engine1600','engine2400','engine3200']:
    x=df[(df.stage=='primary')&(df['filter']==f)&(df.case==case)][m];ys.append(x.mean()*1000);err.append(stats.sem(x)*stats.t.ppf(.975,len(x)-1)*1000)
   ax.errorbar([0,800,1600,2400,3200],ys,yerr=err,label=LABEL[f],color=COLOR[f],marker='o',ms=3,lw=.8)
 axs[0].set_ylabel('Heave RMS [mm]');axs[1].set_ylabel('Horizontal RMS [mm]');axs[1].set_xlabel('Engine speed [rpm; 0 = off]');axs[0].legend(ncol=2);save(fig,out,'engine_degradation')
 d=trace(results,'roundtrip_position');fig,axs=plt.subplots(3,1,figsize=(7.1,3.6),sharex=True)
 axs[0].plot(d.t,d.z_true,color='black',lw=.8,label='True heave')
 for f in LABEL:
  ez=d['z_'+f]-d.z_true;eh=np.hypot(d['x_'+f]-d.x_true,d['y_'+f]-d.y_true)
  axs[1].plot(d.t,ez,color=COLOR[f],lw=.6,label=LABEL[f]);axs[2].plot(d.t,eh,color=COLOR[f],lw=.6)
 for ax in axs:
  for t in [400,800]:ax.axvline(t,color='black',lw=.6,ls=':')
 axs[0].set_ylabel('Heave [m]');axs[1].set_ylabel('Heave error [m]');axs[2].set_ylabel('Horizontal error [m]');axs[2].set_xlabel('Time [s]');axs[1].legend(ncol=2);save(fig,out,'roundtrip')
 fig,axs=plt.subplots(2,1,figsize=(7.1,3.6))
 d=trace(results,'direction_step_position')
 for f in LABEL:axs[0].plot(d.t,np.hypot(d['x_'+f]-d.x_true,d['y_'+f]-d.y_true),lw=.65,color=COLOR[f],label=LABEL[f])
 axs[0].axvline(600,color='black',ls=':',lw=.8);axs[0].set_ylabel('Direction step: H error [m]');axs[0].legend(ncol=2)
 d=trace(results,'mag_outage_attitude')
 for f in LABEL:axs[1].plot(d.t,d['dyaw_'+f],lw=.65,color=COLOR[f])
 axs[1].axvspan(400,700,color='#999999',alpha=.18);axs[1].set_ylabel('Outage: heading error [deg]');axs[1].set_xlabel('Time [s]');save(fig,out,'direction_outage')
def ablations_and_consistency(results,out,df,ab):
 """Plot coupling ablations and uncertainty consistency diagnostics."""
 variants=['TCSM:off','TCSM/real','CSM/real','TCSM/noreset','TCSM:noconf/ext-true','TCSM:noconf/ext-true/real','TCSM/ext-flip','TCSM:noconf/ext-flip']
 short=['Coupling off','Real only (TFG)','Real only (QMEKF)','No reset transport','Known coupling','Known, real only','Wrong + confidence','Wrong, unchecked']
 fig,axs=plt.subplots(2,1,figsize=(3.45,4.7),sharex=True)
 for j,case in enumerate(['long_crested','mixed','direction_step','rapid_rotation']):
  for ax,m in zip(axs,['rms_z','rms_h']):
   s=ab[(ab.case==case)&(ab.metric==m)].set_index('filter').reindex(variants);x=np.arange(len(variants))+(j-1.5)*.15
   ax.errorbar(x,s.percent,yerr=[s.percent-s.ci_low,s.ci_high-s.percent],fmt='o',ms=2.3,lw=.6,label=CASE[case])
 for ax in axs:ax.axhline(0,color='black',lw=.6)
 axs[0].set_ylabel('Heave change [%]');axs[1].set_ylabel('Horizontal change [%]');axs[1].set_xticks(range(len(short)),short,rotation=50,ha='right');axs[0].legend(ncol=2);save(fig,out,'ablations')
 cases=['jonswap','calm','mixed','weak_heave','rapid_rotation','mag_outage','colored_noise'];fig,axs=plt.subplots(2,1,figsize=(3.45,4.2),sharex=True)
 for i,f in enumerate(LABEL):
  x=np.arange(len(cases))+(i-1.5)*.15
  gs=[df[(df.stage=='primary')&(df['filter']==f)&(df.case==case)] for case in cases]
  axs[0].plot(x,[g.nees9_phys.mean() for g in gs],'o',ms=3,color=COLOR[f],label=LABEL[f])
  axs[1].plot(x,[100*g[['cov95_x','cov95_y','cov95_z']].mean().mean() for g in gs],'o',ms=3,color=COLOR[f])
 axs[0].axhline(9,color='black',ls=':',lw=.7);axs[0].set_yscale('log');axs[0].set_ylabel('Joint physical NEES (9 states)');axs[0].legend(ncol=2)
 axs[1].axhline(95,color='black',ls=':',lw=.7);axs[1].set_ylabel('Mean marginal coverage [%]');axs[1].set_xticks(range(len(cases)),[CASE[c] for c in cases],rotation=40,ha='right');save(fig,out,'consistency')
# Reduced-heave diagnosis rows (filter, label); cap0 is the legacy PSD cap.
WEAK=(('A41@C','Q independent'),('T41@C','T independent'),('CSM:off','Q off'),('TCSM:off','T off'),
 ('CSM','Q shared'),('TCSM','T shared'),('CSM/cap0','Q legacy cap'),('TCSM/cap0','T legacy cap'),
 ('TCSM/cap0.7','T cap 0.7'),('TCSM/cap0.9','T cap 0.9'),('CSM/real','Q real only'),('TCSM/real','T real only'),
 ('CSM:noconf','Q unchecked'),('TCSM:noconf','T unchecked'),('CSM:noconf/ext-true','Q known'),
 ('TCSM:noconf/ext-true','T known'),('TCSM/noreset','T no reset'))
def supplemental_figures(results,out):
 """Plot the replay, reset, transition and weak-heave supplements."""
 d=pd.read_csv(results/'weak_runs.csv');names=[f for f,_ in WEAK];labels=[x for _,x in WEAK]
 fig,axs=plt.subplots(2,1,figsize=(7.1,4.0),sharex=True)
 for ax,m,scale,label in zip(axs,['rms_z','yaw_rms_deg'],[1000,1],['Heave RMS [mm]','Heading RMS [deg]']):
  gs=[d[d['filter']==f][m]*scale for f in names]
  ax.bar(np.arange(len(names)),[g.mean() for g in gs],color=['#0072b2' if f.startswith(('CSM','A41')) else '#d55e00' for f in names],alpha=.8)
  ax.errorbar(np.arange(len(names)),[g.mean() for g in gs],yerr=[stats.sem(g)*stats.t.ppf(.975,len(g)-1) for g in gs],fmt='none',color='black',lw=.7,capsize=2)
  ax.set_ylabel(label)
 axs[1].axhline(30,color='black',ls=':',lw=.7);axs[1].set_yscale('log');axs[1].set_xticks(range(len(names)),labels,rotation=45,ha='right');save(fig,out,'weak_diagnostics')
 d=trace(results,'weak_attitude');fig,axs=plt.subplots(2,1,figsize=(7.1,3.5),sharex=True)
 for f in ['CSM','TCSM','CSM/cap0','TCSM/cap0']:
  color=COLOR['TCSM' if f.startswith('TCSM') else 'CSM'];ls='--' if '/cap0' in f else '-'
  label=LABEL['TCSM' if f.startswith('TCSM') else 'CSM']+(', legacy cap' if '/cap0' in f else '')
  axs[0].plot(d.t,d['dyaw_'+f],color=color,ls=ls,lw=.7,label=label)
  c=trace(results,'weak_confidence_'+f.replace('/','_'))
  axs[1].plot(c.t,c.confidence,color=color,ls=ls,lw=.7)
 axs[0].set_ylabel('Heading error [deg]');axs[0].legend(ncol=2);axs[1].set_ylabel('Applied confidence');axs[1].set_xlabel('Time [s]');save(fig,out,'weak_feedback')

# Dose-response series: (filter, label, colour, line, marker); cap0 = legacy cap.
SWEEP=(('TCSM','T shared',COLOR['TCSM'],'-','s'),('CSM','Q shared',COLOR['CSM'],'-','o'),
 ('TCSM/cap0','T legacy cap',COLOR['TCSM'],'--','s'),('CSM/cap0','Q legacy cap',COLOR['CSM'],'--','o'))
def load_heave_sweep(results):
 """Committed sweep scores, protocol, analysis and default-spreading heave factors."""
 runs=pd.read_csv(results/'heave_sweep_runs.csv');p=json.loads((results/'heave_sweep_protocol.json').read_text())
 d=runs[runs.stage=='dose'];iv=runs[runs.stage=='sf2'];ip=p['intervention']
 if len(d)!=len(p['cells'])*len(p['seeds'])*len(p['filters']) or len(iv)!=len(ip['cells'])*len(p['seeds'])*len(ip['filters']):
  raise RuntimeError('Incomplete heave sweep')
 dose={c:float(re.search(r'hscale([0-9.]+)',s)[1]) if 'hscale' in s else 1. for c,s in p['cells'].items() if '+sp' not in s}
 a=json.loads((results/'heave_sweep_analysis.json').read_text())
 return d,iv,p,a,dose
def paired_pct(d,cell,f,b,m):
 """Paired percent change of filter f against b (mean and 95% half-width)."""
 v=d[d.cell==cell].pivot(index='seed',columns='filter',values=m);x=100*(v[f]/v[b]-1)
 return x.mean(),stats.sem(x)*stats.t.ppf(.975,len(x)-1)
def heave_sweep(results,out):
 """Plot heave, horizontal and heading against the heave RAO for both cap policies."""
 d,_,p,_,dose=load_heave_sweep(results);ks=list(dose.values())
 fig,axs=plt.subplots(3,1,figsize=(3.45,5.2),sharex=True)
 for f,label,color,ls,marker in SWEEP:
  for ax,m in zip(axs[:2],('rms_z','rms_h')):
   r=[paired_pct(d,c,f,'T41@C',m) for c in dose]
   ax.errorbar(ks,[x for x,_ in r],yerr=[h for _,h in r],color=color,ls=ls,marker=marker,ms=3,lw=.8,capsize=1.5,label=label)
  q=np.array([np.percentile(d[(d.cell==c)&(d['filter']==f)].yaw_rms_deg,[25,50,75]) for c in dose])
  axs[2].errorbar(ks,q[:,1],yerr=[q[:,1]-q[:,0],q[:,2]-q[:,1]],color=color,ls=ls,marker=marker,ms=3,lw=.8,capsize=1.5)
 q=[d[(d.cell==c)&(d['filter']=='T41@C')].yaw_rms_deg.median() for c in dose]
 axs[2].plot(ks,q,color=COLOR['T41@C'],marker='^',ms=3,lw=.8,label='T independent')
 for ax in axs[:2]:ax.axhline(0,color='black',lw=.6);ax.set_yscale('symlog',linthresh=20)
 axs[0].set_ylabel('Heave vs T indep. [%]');axs[1].set_ylabel('Horizontal vs T indep. [%]');axs[0].legend(fontsize=5.5,ncol=2,loc='upper left')
 axs[2].axhline(p['heading_failure_deg'],color='black',ls=':',lw=.7);axs[2].set_yscale('log');axs[2].set_ylabel('Heading RMS [deg]');axs[2].legend(fontsize=5.5,loc='upper left')
 axs[2].set_xscale('log');axs[2].set_xlim(1.12,.135);axs[2].set_xticks(ks,[f'{x:g}' for x in ks]);axs[2].minorticks_off()
 axs[2].set_xlabel('Heave RAO factor $k$ (horizontal motion unchanged)');save(fig,out,'heave_sweep')

def heave_sweep_claims(results):
 """Tabulate the sweep and derive its discussion only from measured cells."""
 d,iv,p,a,dose=load_heave_sweep(results);n=len(p['seeds']);fail=p['heading_failure_deg'];cells=list(p['cells'])
 cap=a['cap_binding'];harm=a['legacy_heave_vs_real_only']
 label={c:f'{x:g}' for c,x in dose.items()}
 label.update({c:'0.15 broad' for c in cells if c.endswith('broad')});label.update({c:'0.15 long' for c in cells if c.endswith('long')})
 mean=lambda c,f,m:d[(d.cell==c)&(d['filter']==f)][m].mean()
 coh=lambda c,f:(d[(d.cell==c)&(d['filter']==f)].cpl_coh_xz+d[(d.cell==c)&(d['filter']==f)].cpl_coh_yz).mean()
 rows=[r'\begin{table*}[!t]\centering',
  r'\caption{Heave-RAO dose response, '+str(n)+r' paired seeds per cell. $k$ scales the heave RAO only; horizontal motion and every sensor stream are unchanged within a seed; broad and long-crested cells change only the $\cos^{2s}$ spreading ($s=1$, $s=10^4$; default $s=4$). $\gamma^2$: applied multiple coherence of the known diagnostic prior with the legacy cap (0.999 ceiling) and the default cap ('+f"{p['fixed_policy']['cap_coherence']:g}"+r'). RMS errors are means; heading is the median RMS with runs above '+f'{fail:g}'+r'$^\circ$ in parentheses. Legacy: cross term scaled down at the cap.}',
  r'\footnotesize\begin{tabular}{lrrrrrrrrrrrr}\toprule',
  r'&\multicolumn{2}{c}{$\gamma^2$}&\multicolumn{5}{c}{Heave RMS [mm]}&\multicolumn{3}{c}{Horizontal RMS [mm]}&\multicolumn{2}{c}{Heading [deg]}\\',
  r'$k$&legacy&default&T ind.&Q legacy&T legacy&Q shared&T shared&T ind.&Q shared&T shared&T legacy&T shared\\\midrule']
 for c in cells:
  vals=[f"{coh(c,'CSM:noconf/ext-true/cap0'):.3f}",f"{coh(c,'CSM:noconf/ext-true'):.3f}"]
  vals+=[fmt(1000*mean(c,f,'rms_z'),1) for f in ('T41@C','CSM/cap0','TCSM/cap0','CSM','TCSM')]
  vals+=[fmt(1000*mean(c,f,'rms_h'),1) for f in ('T41@C','CSM','TCSM')]
  for f in ('TCSM/cap0','TCSM'):
   x=d[(d.cell==c)&(d['filter']==f)].yaw_rms_deg;k=int((x>fail).sum());vals.append(f'{x.median():.1f}'+(f' ({k})' if k else ''))
  rows.append(label[c]+'&'+'&'.join(vals)+r'\\')
 rows.append(r'\bottomrule\end{tabular}\label{tab:heave-sweep}\end{table*}')
 iv_a=a['intervention_sf2'];sf=p['intervention']['s_factor'];iv_cells=p['intervention']['cells']
 rows.extend([r'\begin{table}[!t]\centering',
  r'\caption{Legacy-cap intervention: the cells of Table~\ref{tab:heave-sweep} near the cap with only the horizontal prior factor raised to $S='+f'{sf:g}'+r'$ (same seeds). Heave change relative to legacy real-only coupling at the same $S$, with paired 95\% intervals; mean heading RMS with runs above '+f'{fail:g}'+r'$^\circ$ in parentheses.}',
  r'\footnotesize\begin{tabular}{lrllrr}\toprule $k$&$\gamma^2$&Q legacy [\%]&T legacy [\%]&Q [deg]&T [deg]\\\midrule'])
 for c in iv_cells:
  e=iv_a[c];vals=[f"{e[f]['heave_vs_real_pct']:+.1f} [{e[f]['ci_low']:+.1f}, {e[f]['ci_high']:+.1f}]" for f in ('CSM/cap0/sf2','TCSM/cap0/sf2')]
  vals+=[f"{e[f]['heading_mean_deg']:.1f}"+(f" ({e[f]['heading_over30']})" if e[f]['heading_over30'] else '') for f in ('CSM/cap0/sf2','TCSM/cap0/sf2')]
  rows.append(f"{dose[c]:g}&{e['known_coherence']:.3f}&"+'&'.join(vals)+r'\\')
 rows.append(r'\bottomrule\end{tabular}\label{tab:heave-sf2}\end{table}')
 rep=a['reproduction'];rows_rep=sum(v['rows'] for v in rep.values());exact=all(v['max_relative_difference']==0 for v in rep.values())
 capped=[c for c in dose if cap[c]['known_at_cap']];first_cap=max(dose[c] for c in capped) if capped else None
 uncapped=[c for c in dose if not cap[c]['known_at_cap']]
 hq=harm['CSM/cap0']['first_harm_cell'];onset=a['heading_onset']
 kcap=(cap['h1.00']['known_legacy']/(sf*sf*0.999))**.5
 below=[c for c in iv_cells if iv_a[c]['known_coherence']<0.999-5e-4];at=[c for c in iv_cells if c not in below]
 shift=all(iv_a[c][f]['ci_high']<0 for c in below for f in ('CSM/cap0/sf2','TCSM/cap0/sf2')) and all(iv_a[c][f]['ci_low']>0 for c in at for f in ('CSM/cap0/sf2','TCSM/cap0/sf2'))
 four=('A41@C','T41@C','CSM','TCSM')
 excess=lambda cs,m:max(0.,max(100*(mean(c,'TCSM',m)/min(mean(c,f,m) for f in four)-1) for c in cs))
 margin=lambda best,m:'' if best==len(cells) else f" (at most {excess(cells,m):.2f}\\% above the lowest elsewhere)"
 best_z=sum(min(('A41@C','T41@C','CSM','TCSM'),key=lambda f:mean(c,f,'rms_z'))=='TCSM' for c in cells)
 best_h=sum(min(('A41@C','T41@C','CSM','TCSM'),key=lambda f:mean(c,f,'rms_h'))=='TCSM' for c in cells)
 rz=[paired_pct(d,c,'TCSM','T41@C','rms_z')[0] for c in dose];rh=[paired_pct(d,c,'TCSM','T41@C','rms_h')[0] for c in dose]
 fixed_fail=int(sum((d[(d['filter'].isin(['CSM','TCSM']))].yaw_rms_deg>fail)))
 pct=lambda r:f"{r['percent']:+.1f}\\% [{r['ci_low']:+.1f}, {r['ci_high']:+.1f}]"
 legacy_coh=', '.join(f"{cap[c]['known_legacy']:.3f}" for c in uncapped[1:]);legacy_k=', '.join(label[c] for c in uncapped[1:])
 text=(f"A heave-RAO dose response (Tables~\\ref{{tab:heave-sweep}} and~\\ref{{tab:heave-sf2}}, Fig.~\\ref{{fig:heave-sweep}}) repeats the first {n} primary seeds with the heave RAO scaled by $k$; horizontal motion and every sensor stream are unchanged within a seed. "
  +(f"Its repeated $k=1$ and $k=0.15$ cells reproduce the heave, horizontal, heading, tilt and NEES scores of all {rows_rep} overlapping primary and diagnosis runs exactly. " if exact else "")
  +f"With the legacy cap the known prior's applied multiple coherence is {cap['h1.00']['known_legacy']:.3f} at $k=1$ and grows as $1/k^2$ "
  +f"({legacy_coh} at $k={legacy_k}$), "
  +(f"reaching the 0.999 cap for $k\\le{first_cap:g}$. " if capped else "without reaching the cap. ")
  +f"The legacy shared-QMEKF heave change relative to legacy real-only coupling is {pct(harm['CSM/cap0']['h1.00'])} at $k=1$ and {pct(harm['CSM/cap0'][uncapped[-1]])} at $k={label[uncapped[-1]]}$, "
  +(f"but {pct(harm['CSM/cap0'][hq])} at $k={dose[hq]:g}$, where the cap binds; " if hq else "and never turns into harm; ")
  +(f"the legacy heading interval first exceeds the independent filter by one degree at $k={dose[onset['CSM/cap0']]:g}$ (QMEKF) and $k={dose[onset['TCSM/cap0']]:g}$ (TFG). " if onset['CSM/cap0'] and onset['TCSM/cap0'] else "")
  +(f"Raising only the horizontal prior factor to {sf:g} moves the predicted binding boundary to $k={kcap:.2f}$ and the observed one to between {min(dose[c] for c in below):g} and {max(dose[c] for c in at):g}. " if shift and at else "The intervention does not move the boundary as predicted. ")
  +"The capped legacy prior is nearly singular along the coupling direction: 99.9\\% of the along-wave horizontal variance is tied to heave through an understated regression, while the cross-wave direction keeps its full variance. "
  +f"The default policy keeps the regression and holds the applied coherence at {p['fixed_policy']['cap_coherence']:g} by enlarging the horizontal prior along the coupling direction only. "
  +f"Across the {len(cells)} cells, shared TFG then has the lowest mean heave RMS of the four filters in {best_z}{margin(best_z,'rms_z')} and the lowest horizontal RMS in {best_h}{margin(best_h,'rms_h')}; "
  +f"relative to independent TFG its heave changes by {min(rz):+.1f}\\% to {max(rz):+.1f}\\% and its horizontal error by {min(rh):+.1f}\\% to {max(rh):+.1f}\\% over the default-spreading factors, "
  +f"and {fixed_fail} shared runs exceed {fail:g} degrees heading.")
 return rows,text,dict(cap_shift=shift,best_z=best_z,best_h=best_h,cells=len(cells))

def supplemental_tables(results,main):
 """Summarize supplementary scores without pooling independent wave records."""
 weak=pd.read_csv(results/'weak_summary.csv').set_index('filter');wr=pd.read_csv(results/'weak_runs.csv')
 ds=pd.read_csv(results/'dataset_runs.csv');dp=pd.read_csv(results/'dataset_paired.csv')
 reset=pd.read_csv(results/'reset_comparison.csv');tr=pd.read_csv(results/'transient_runs.csv');tp=pd.read_csv(results/'transient_paired.csv')
 if len(wr)!=12*len(WEAK) or len(ds)!=8*8*4 or len(tr)!=24*3*4 or len(reset)!=6*2*7:
  raise RuntimeError('Incomplete supplementary experiments')
 failures=main.groupby(['case','filter']).agg(n=('seed','size'),heading_over30=('yaw_rms_deg',lambda x:int((x>30).sum())),heading_mean_deg=('yaw_rms_deg','mean'),benchmark_flags=('failed','sum'))
 failures.to_csv(results/'heading_failures.csv',float_format='%.9g')
 rows=[r'\begin{table*}[!t]\centering',r'\caption{Exploratory reduced-heave diagnosis, 12 paired seeds. Q = QMEKF, T = TFG. Legacy cap: the cross term is scaled down at the PSD cap; cap 0.7/0.9 bracket the default target 0.8. Known coupling is a diagnostic input. The 30-degree count is separate from the inherited divergence flag.}',r'\footnotesize\begin{tabular}{lrrrrrrr}\toprule Variant&Heave [mm]&H [mm]&Heading [deg]&Heading $>30^\circ$&NEES$_9$&Confidence&Shrinks\\\midrule']
 for f,label in WEAK:
  r=weak.loc[f]
  vals=[fmt(r.heave_m*1000),fmt(r.horizontal_m*1000),fmt(r.heading_deg),f'{int(r.heading_over30)}/12',fmt(r.nees9,0),fmt(r.mean_confidence),fmt(r.shrinks,0)]
  rows.append(tex(label)+'&'+'&'.join(vals)+r'\\')
 rows.extend([r'\bottomrule\end{tabular}\label{tab:weak}\end{table*}'])
 rows.extend([r'\begin{table*}[!t]\centering',r'\caption{Pinned vessel-RAO replay, eight sensor/initialization seeds conditional on each fixed wave record. Heave RMS in cm. T/Q shared change has a paired 95\% interval; wave realizations are not independently resampled.}',r'\footnotesize\begin{tabular}{lrrrrl}\toprule Record&Q ind.&T ind.&Q shared&T shared&T/Q shared [\%]\\\midrule'])
 for record,g in ds.groupby('record',sort=True):
  typ='JONSWAP' if 'jonswap' in record else 'PM-Stokes';height=re.search(r'_H([0-9.]+)_L',record)[1]
  vals=[100*g[g.method==f].disp_z_rms_m.mean() for f in ['Q independent','T independent','Q shared','T shared']]
  r=dp[(dp.record==record)&(dp.metric=='disp_z_rms_m')].iloc[0]
  rows.append(f'{typ}, H={height} m&'+'&'.join(fmt(x,3) for x in vals)+f'&{r.percent:+.4f} [{r.ci_low:+.4f}, {r.ci_high:+.4f}]'+r'\\')
 rows.extend([r'\bottomrule\end{tabular}\end{table*}'])
 rows.extend([r'\begin{table*}[!t]\centering',r'\caption{Reset diagnostic: change when covariance transport is omitted, percent with paired 95\% intervals (12 seeds). The default retains full transport.}',r'\footnotesize\begin{tabular}{llrr}\toprule Condition&Prior&Heave [\%]&Horizontal [\%]\\\midrule'])
 for (case,geometry),g in reset.groupby(['case','geometry']):
  vals=[]
  for m in ['rms_z','rms_h']:
   r=g[g.metric==m].iloc[0];vals.append(f'{r.no_transport_pct:+.3f} [{r.ci_low:+.3f}, {r.ci_high:+.3f}]')
  rows.append(tex(CASE[case])+'&'+geometry+'&'+'&'.join(vals)+r'\\')
 rows.extend([r'\bottomrule\end{tabular}\end{table*}'])
 rows.extend([r'\begin{table*}[!t]\centering',r'\caption{Transition recovery (24 seeds). Mean times in seconds; paired TFG/shared minus QMEKF/shared intervals. Recovery is relative to each filter\textquotesingle s own post-event RMS and does not establish absolute accuracy.}',r'\footnotesize\begin{tabular}{lrrrrll}\toprule Event&Q heave&T heave&Q horizontal&T horizontal&Heave difference [95\% CI]&H difference [95\% CI]\\\midrule'])
 for (case,event),g in tr.groupby(['case','event']):
  vals=[]
  for m in ['recovery_z_s','recovery_h_s']:
   vals.extend([fmt(g[g['filter']==f][m].mean()) for f in ['CSM','TCSM']])
  intervals=[]
  for m in ['recovery_z_s','recovery_h_s']:
   r=tp[(tp.case==case)&(tp.event==event)&(tp['filter']=='TCSM')&(tp.metric==m)].iloc[0]
   intervals.append(f'{r.delta:+.2f} [{r.ci_low:+.2f}, {r.ci_high:+.2f}]')
  rows.append(tex(CASE[case]+' '+event)+'&'+'&'.join(vals+intervals)+r'\\')
 rows.extend([r'\bottomrule\end{tabular}\end{table*}'])
 rows.extend([r'\begin{table}[!t]\centering',r'\caption{Nominal moving startup and magnetic-outage recovery. Times omit censored runs; counts are shown explicitly.}',r'\footnotesize\begin{tabular}{lrrrr}\toprule &Q ind.&T ind.&Q shared&T shared\\\midrule'])
 for label,case,col,scale in [('Startup heave [cm]','jonswap','su_rms_z',100),('Startup tilt [deg]','jonswap','su_tilt_rms_deg',1),('Startup yaw [deg]','jonswap','su_yaw_rms_deg',1),('Heave convergence [s]','jonswap','conv_z',1),('Outage yaw [deg]','mag_outage','yaw_gap_rms_deg',1),('Post-outage yaw [deg]','mag_outage','yaw_post_rms_deg',1),('Yaw recovery [s]','mag_outage','yaw_recover_s',1)]:
  gs=[main[(main.case==case)&(main['filter']==f)][col] for f in LABEL]
  rows.append(label+'&'+'&'.join(fmt(g.mean()*scale) for g in gs)+r'\\')
  if col in ['conv_z','yaw_recover_s']:rows.append('Observed / 24'+'&'+'&'.join(str(int(g.notna().sum())) for g in gs)+r'\\')
 rows.extend([r'\bottomrule\end{tabular}\end{table}'])
 cens=int(tr[['censored_z','censored_h']].to_numpy().sum())
 diagnostic={'weak':weak,'heading_failures':failures,'transient_censored':cens,'replay_heave_min':dp[dp.metric=='disp_z_rms_m'].percent.min(),'replay_heave_max':dp[dp.metric=='disp_z_rms_m'].percent.max()}
 return rows,diagnostic
def tables_and_macros(results,df,pairs,agg,ab):
 """Validate the complete study and derive article tables, claims and input hashes."""
 protocol=json.loads((results/'protocol.json').read_text());ns=len(protocol['primary_seeds']);na=len(protocol['ablation_seeds']);nc=len(protocol['scenarios'])
 complete=all(len(df[(df.stage=='primary')&(df.case==c)])==ns*4 for c in protocol['scenarios']) and len(df[df.stage=='ablation'])==na*len(protocol['ablation_scenarios'])*len(protocol['ablation_filters'])
 if not complete:raise RuntimeError('The scientific article requires the complete registered experiment.')
 target=pd.read_csv(results/'esp32s3_size_stack.csv').set_index('variant')
 def aggregate(m):
  """Select the condition-balanced shared-TFG geometry comparison for one metric."""
  return agg[(agg.metric==m)&(agg.comparison=='TFG geometry, shared')].iloc[0]
 z=aggregate('rms_z');h=aggregate('rms_h');bad=int(df.failed.sum());main=df[df.stage=='primary']
 extra_rows,diagnostic=supplemental_tables(results,main);sweep_rows,sweep_text,sweep=heave_sweep_claims(results)
 htime=pd.read_csv(results/'host_timing.csv');jtime=pd.read_csv(results/'host_jitter.csv')
 # Validate the fixed-model timing input even though this table uses commits.
 pd.read_csv(results/'host_fixed.csv')
 t=htime.groupby('core').median(numeric_only=True);jt=jtime.groupby('core').median(numeric_only=True)
 def eff(r):
  """Render a measured aggregate effect and its paired confidence interval."""
  return f"{r.geomean_pct:+.2f}\\% (95\\% CI {r.ci_low_pct:+.2f} to {r.ci_high_pct:+.2f}\\%)"
 best=pairs[(pairs.comparison=='TFG geometry, shared')&pairs.metric.isin(['rms_z','rms_h'])]
 practical=best[(best.pct_ci_high<-5)&(best.p_holm<.05)]
 adverse=best[(best.percent>0)&(best.p_holm<.05)]
 adverse_txt=f", at most {adverse.percent.max():+.2f}\\%" if len(adverse) else ''
 macro={
 'ResultsFigureDirectory':os.path.relpath(results/'figures',DOC).replace(os.sep,'/'),
 'PrimarySeeds':str(ns),'AblationSeeds':str(na),'ScenarioCount':str(nc),
 'AbstractResults':f"The condition-balanced TFG change relative to shared QMEKF is {eff(z)} for heave and {eff(h)} for horizontal position. There are {bad} flagged estimator trajectories across all registered operational and diagnostic cells.",
 'RegressionResults':"The accompanying regression record lists the original and post-refactor gates and exact paired numerical differences. Existing thresholds are retained. The full-wrapper coupling-off smoke uses the existing TFG/CARMA 1\\% agreement criterion; the process-law identities use $10^{-11}$ relative tolerance.",
 'ResultsDiscussion':f"Condition-balanced heave changes by {eff(z)} and horizontal position by {eff(h)} when replacing shared QMEKF with shared TFG. {len(practical) or 'No'} scenario/endpoint{'' if len(practical)==1 else 's'} meet{'s' if len(practical)==1 else ''} both the predeclared 5\\% improvement margin and the Holm significance criterion; {len(adverse)} show{'s' if len(adverse)==1 else ''} statistically supported worsening{adverse_txt}. The complete per-condition table prevents small aggregate changes from hiding different behavior in weak excitation, maneuvers or outages. The interaction table and figure evaluate nonadditivity directly. There are {bad} flagged trajectories; their identities and the remaining consistency defects are retained below. Exact conditional-prior algebra does not make the inherited wave/bias/attitude approximation statistically calibrated.",
 'PerformanceResults':f"On this host, the median of three repetitions is {t.loc['TCSM','total_us']:.2f} us per 200 Hz core step for shared TFG, versus {t.loc['CSM','total_us']:.2f} us for shared QMEKF (ratio {t.loc['TCSM','total_us']/t.loc['CSM','total_us']:.2f}). With 25\\% interval jitter the corresponding times are {jt.loc['TCSM','total_us']:.2f} and {jt.loc['CSM','total_us']:.2f} us. Timings include 0.1 s model commits and magnetic updates at 25 Hz; percentiles and maxima are reported separately, and host scheduling can affect their tails.",
 'FirmwareResults':f"The Xtensa probe gives {int(target.loc['TCSM','sizeof_qmekf'])} bytes for the shared TFG core versus {int(target.loc['CSM','sizeof_qmekf'])} for shared QMEKF. Their largest emitted function stack frames are {int(target.loc['TCSM','max_stack_frame'])} and {int(target.loc['CSM','max_stack_frame'])} bytes; the static update-chain lower bounds are {int(target.loc['TCSM','update_chain_stack'])} and {int(target.loc['CSM','update_chain_stack'])} bytes. The compiler used for this probe is Espressif GCC esp-14.2.0\\_20241119 with \\texttt{{-Os}}. A full Arduino firmware build is recorded separately in the reproducibility bundle.",
 }
 verification=(results/'mathematical_verification.txt').read_text()
 def largest_check(fragment):
  """Require and summarize a named mathematical verification diagnostic."""
  vals=[float(x) for name,x in re.findall(r'CHECK (.*?)value=([0-9.eE+-]+)',verification) if fragment in name]
  if not vals:raise RuntimeError('Missing mathematical verification: '+fragment)
  return max(vals)
 def sci(x):
  """Render a measured numerical discrepancy in scientific LaTeX notation."""
  a,b=f'{x:.2e}'.split('e');return a+r'\times10^{'+str(int(b))+'}'
 pred=largest_check('prediction Jacobian');cov=largest_check('Hermitian Qd');coords=largest_check('adaptive coordinate covariance')
 macro['VerificationResults']=f"Measured maximum prediction finite-difference discrepancy is ${sci(pred)}$; the relative independent continuous-covariance discrepancy is ${sci(cov)}$. The anisotropic model-coordinate covariance discrepancy is ${sci(coords)}$. Allocation-guarded float execution explicitly forces a falling-confidence fallback and remains finite with symmetric covariance; round-off eigenvalues are recorded rather than hidden."
 nominal=main[main.case=='jonswap'].groupby('filter').mean(numeric_only=True)
 macro['ResultsDiscussion'] += f" In the nominal case, joint physical NEES is {nominal.loc['TCSM','nees9_phys']:.1f} for shared TFG and {nominal.loc['CSM','nees9_phys']:.1f} for shared QMEKF, against a 9-state reference mean of 9. This diagnostic is not repaired by the geometry change."
 weak=main[main.case=='weak_heave'].groupby('filter').mean(numeric_only=True)
 failures=diagnostic['heading_failures'];w=diagnostic['weak'];nw=int(w.loc['CSM','n'])
 over=lambda f:int(failures.loc[('weak_heave',f),'heading_over30'])
 flipped=ab[(ab.case=='rapid_rotation')&(ab['filter']=='TCSM:noconf/ext-flip')&(ab.metric=='rms_h')]
 if len(flipped)!=1:raise RuntimeError('Missing or ambiguous reversed-coupling ablation')
 flipped_change=float(flipped.iloc[0].percent)
 def comparison(m,label):
  """Condition-balanced aggregate row for one metric and comparison."""
  return agg[(agg.metric==m)&(agg.comparison==label)].iloc[0]
 prior_z,prior_h=comparison('rms_z','Shared prior, TFG'),comparison('rms_h','Shared prior, TFG')
 best_z,best_h=comparison('rms_z','TFG shared vs QMEKF independent'),comparison('rms_h','TFG shared vs QMEKF independent')
 mm=lambda f,m:1000*weak.loc[f,m]
 macro['AbstractResults']+=(f" Relative to independent-axis TFG, the shared prior changes heave by {eff(prior_z)} and horizontal position by {eff(prior_h)}; "
  f"relative to independent QMEKF, shared TFG changes them by {eff(best_z)} and {eff(best_h)}. "
  f"Under reduced heave shared TFG has {mm('TCSM','rms_z'):.1f}/{mm('TCSM','rms_h'):.1f} mm heave/horizontal error against {mm('T41@C','rms_z'):.1f}/{mm('T41@C','rms_h'):.1f} mm independently, "
  f"with {over('TCSM')}/{ns} heading failures; the legacy excitation cap gives {1000*w.loc['TCSM/cap0','heave_m']:.1f}/{1000*w.loc['TCSM/cap0','horizontal_m']:.1f} mm and {int(w.loc['TCSM/cap0','heading_over30'])}/{nw} heading failures in a paired {nw}-seed diagnosis.")
 real_worse=min(w.loc['CSM/real','heave_m'],w.loc['TCSM/real','heave_m'])>max(w.loc['A41@C','heave_m'],w.loc['T41@C','heave_m'])
 real_vs_independent=', worse than the independent filters, so the quadrature term is required here' if real_worse else ''
 macro['FailureDiscussion']=(f"Reduced heave is the condition that separates the cap policies. With the default policy, mean heave/horizontal errors are "
  f"{mm('CSM','rms_z'):.1f}/{mm('CSM','rms_h'):.1f} mm for shared QMEKF and {mm('TCSM','rms_z'):.1f}/{mm('TCSM','rms_h'):.1f} mm for shared TFG, against "
  f"{mm('A41@C','rms_z'):.1f}/{mm('A41@C','rms_h'):.1f} and {mm('T41@C','rms_z'):.1f}/{mm('T41@C','rms_h'):.1f} mm independently; "
  f"{over('CSM')} and {over('TCSM')} of {ns} shared runs exceed 30 degrees RMS heading. "
  f"In the paired {nw}-seed diagnosis (Table~\\ref{{tab:weak}}) the legacy cap, which scales the cross term down, gives {1000*w.loc['CSM/cap0','heave_m']:.1f} and {1000*w.loc['TCSM/cap0','heave_m']:.1f} mm heave, "
  f"{w.loc['CSM/cap0','heading_deg']:.1f} and {w.loc['TCSM/cap0','heading_deg']:.1f} degrees mean heading RMS and {int(w.loc['CSM/cap0','heading_over30'])} and {int(w.loc['TCSM/cap0','heading_over30'])} heading failures. "
  f"Cap targets 0.7 and 0.9 give {1000*w.loc['TCSM/cap0.7','heave_m']:.1f} and {1000*w.loc['TCSM/cap0.9','heave_m']:.1f} mm heave for shared TFG, against {1000*w.loc['TCSM','heave_m']:.1f} mm at the default 0.8. "
  f"Without confidence weighting the estimated coupling gives {1000*w.loc['CSM:noconf','heave_m']:.1f} and {1000*w.loc['TCSM:noconf','heave_m']:.1f} mm and the known coupling {1000*w.loc['CSM:noconf/ext-true','heave_m']:.1f} and {1000*w.loc['TCSM:noconf/ext-true','heave_m']:.1f} mm, "
  f"so the estimator recovers the known-coupling benefit and confidence weighting gives up part of it in this condition. "
  f"The real-only control maps only the co-spectrum and gives {1000*w.loc['CSM/real','heave_m']:.1f} and {1000*w.loc['TCSM/real','heave_m']:.1f} mm{real_vs_independent}. "
  f"Omitting reset transport gives {1000*w.loc['TCSM/noreset','heave_m']:.1f} mm shared-TFG heave, so the transport is immaterial in this condition.")
 macro['HeaveSweepDiscussion']=sweep_text
 macro['SupplementDiscussion']=f"The pinned eight-record replay gives heave changes between {diagnostic['replay_heave_min']:+.4f}\\% and {diagnostic['replay_heave_max']:+.4f}\\%, too small to justify a geometry switch. These intervals condition on each wave record and resample only sensors/initialization. The 24-seed full-rate transition supplement finds no practically meaningful recovery advantage; {diagnostic['transient_censored']} transition recoveries are censored. Its exported positions have 10-micrometre resolution, so sub-resolution differences are not evidence for useful improvement. The six-condition reset supplement shows small ordinary-case position changes for both geometries. Reversed unchecked coupling degrades horizontal error by {flipped_change:.2f}\\% in rapid rotations, while confidence reduces that diagnostic to the uncoupled response. Confidence is useful but does not guarantee robustness."
 engine=main[main.case=='engine800'].groupby('filter').mean(numeric_only=True)
 turn=main[main.case=='maneuver'].groupby('filter').mean(numeric_only=True)
 colored=main[main.case=='colored_noise'].groupby('filter').mean(numeric_only=True)
 yaw=pairs[(pairs.case=='maneuver')&(pairs.comparison=='TFG geometry, shared')&(pairs.metric=='yaw_rms_deg')].iloc[0]
 def worded(frame,m,what,unit):
  """Shared QMEKF to shared TFG, worded by the sign of the measured change."""
  a,b=frame.loc['CSM',m],frame.loc['TCSM',m]
  return f"{'reduces' if b<a else 'increases'} {what} from {a:.2f} to {b:.2f} {unit}"
 horizontal=pairs[(pairs.case=='colored_noise')&(pairs.comparison=='TFG geometry, shared')&(pairs.metric=='rms_h')].iloc[0]
 bg=pairs[(pairs.case=='maneuver')&(pairs.comparison=='TFG geometry, shared')&(pairs.metric=='bg_rms_dph')].iloc[0]
 macro['SecondaryDiscussion']=(f"Secondary attitude and bias results are condition dependent. At 800 rpm, shared TFG {worded(engine,'tilt_rms_deg','mean tilt RMS','degrees')}, "
  f"{worded(engine,'yaw_rms_deg','heading','degrees')} and {worded(engine,'ba_rms_mg','accelerometer-bias error','mg')}. "
  f"Under colored noise it {worded(colored,'yaw_rms_deg','heading','degrees')}, while horizontal position changes by {horizontal.percent:+.2f}\\% [{horizontal.pct_ci_low:+.2f}, {horizontal.pct_ci_high:+.2f}]. "
  f"During maneuvers it {worded(turn,'yaw_rms_deg','heading','degrees')}: paired absolute change {yaw.delta:+.2f} degrees, marginal 95\\% CI [{yaw.ci_low:+.2f}, {yaw.ci_high:+.2f}]; "
  f"gyro-bias RMS changes by {bg.percent:+.2f}\\% [{bg.pct_ci_low:+.2f}, {bg.pct_ci_high:+.2f}]. "
  "These secondary intervals are exploratory and are not part of the Holm-controlled primary position family.")
 geo_z,geo_h=comparison('rms_z','TFG geometry, shared'),comparison('rms_h','TFG geometry, shared')
 vs_ind=pairs[(pairs.comparison=='TFG shared vs QMEKF independent')&pairs.metric.isin(['rms_z','rms_h'])&(pairs.case!='calm')]
 worst=vs_ind.loc[vs_ind.percent.idxmax()]
 best_all=all(r.ci_high_pct<0 for r in (geo_z,geo_h,prior_z,prior_h,best_z,best_h))
 ratio=t.loc['TCSM','total_us']/t.loc['CSM','total_us']
 macro['Recommendation']=(("Shared TFG is the most accurate of the four estimators in the condition-balanced position comparison: it improves on shared QMEKF, independent TFG and independent QMEKF with paired intervals excluding zero for heave and horizontal position. "
   if best_all else "Shared TFG is not uniformly the most accurate of the four estimators: at least one condition-balanced position comparison has an interval that includes zero or favours another filter. ")
  +f"Most of the gain comes from the shared prior ({eff(prior_z)} heave, {eff(prior_h)} horizontal); the geometry adds {eff(geo_z)} and {eff(geo_h)} at {ratio:.1f} times the host cost of shared QMEKF. "
  +"The shared prior is worth keeping in both geometries with the cap policy; TFG is the better choice only where its small accuracy and attitude gains justify the added computation. "
  +f"The largest condition-specific cost against independent QMEKF is {worst.percent:+.2f}\\% {('heave' if worst.metric=='rms_z' else 'horizontal')} RMS ({tex(CASE[worst.case])}). Nominal joint NEES remains poorly calibrated. Actual ESP32 timing and stack high-water measurements are still required. The shipping OU-III default remains unchanged.")
 firmware=json.loads((results/'firmware_builds.json').read_text())
 builds={r['variant']:r for r in firmware['builds']}
 if any(v not in builds or builds[v]['exit_code'] for v in ['qmekf','tfg','probe']):
  raise RuntimeError('Final article requires successful matched firmware and runtime-probe builds')
 macro['FirmwareResults']+=f" With matched Arduino optimization flags, the complete AtomS3R shared-mode application links at {builds['tfg']['flash_bytes']:,} bytes for TFG and {builds['qmekf']['flash_bytes']:,} bytes for QMEKF; both use {builds['tfg']['static_ram_bytes']:,} bytes of globals. These figures exclude the estimator allocated at construction, whose core size is tabulated separately. The board package is ESP32 3.3.7 with M5Unified 0.2.13 and Espressif GCC 14.2.0\\_20251107. The 200 Hz measurement harness also compiles, but its serial timing/heap/stack results require a connected board."
 DOC.mkdir(parents=True,exist_ok=True);(DOC/'tfg-shared-results.tex-part').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macro.items())+'\n')
 rows=[r'\begin{table*}[!t]\centering',r'\caption{Full-rate primary position errors and paired geometry changes. RMS is in cm; changes are percent, with marginal paired 95\% intervals. Q = QMEKF, T = TFG.}',r'\footnotesize\begin{tabular}{lrrrrll}\toprule',r'&Q independent& T independent& Q shared& T shared&\multicolumn{2}{c}{T shared versus Q shared, \%}\\',r'Condition&\multicolumn{4}{c}{Heave RMS [cm]}&Heave [95\% CI]&Horizontal [95\% CI]\\\midrule']
 for case in CASE:
  g=main[main.case==case].set_index('filter');val=[g.loc[f,'rms_z'].mean()*100 for f in LABEL]
  pz=pairs[(pairs.case==case)&(pairs.metric=='rms_z')&(pairs.comparison=='TFG geometry, shared')].iloc[0];ph=pairs[(pairs.case==case)&(pairs.metric=='rms_h')&(pairs.comparison=='TFG geometry, shared')].iloc[0]
  rows.append(tex(CASE[case])+'&'+'&'.join(fmt(v) for v in val)+'&'+interval(pz)+'&'+interval(ph)+r'\\')
 rows.extend([r'\bottomrule\end{tabular}',r'\label{tab:primary}\end{table*}'])
 rows.extend([r'\begin{table}[!t]\centering\caption{Condition-balanced paired TFG/shared versus QMEKF/shared changes. Seeds are independent; each condition has equal log-ratio weight.}',r'\footnotesize\begin{tabular}{lrr}\toprule Metric&Change [\%]&95\% CI [\%]\\\midrule'])
 for m,label in [('rms_z','Heave'),('rms_h','Horizontal'),('rms_3d','3-D position'),('tilt_rms_deg','Tilt'),('yaw_rms_deg','Heading'),('bg_rms_dph','Gyro bias'),('ba_rms_mg','Accel. bias')]:
  r=aggregate(m);rows.append(f'{label}&{r.geomean_pct:+.3f}&[{r.ci_low_pct:+.3f}, {r.ci_high_pct:+.3f}]'+r'\\')
 rows.extend([r'\bottomrule\end{tabular}\end{table}'])
 rows.extend([r'\begin{table}[!t]\centering\caption{Host core and Xtensa compile-time measurements. Static stack bounds exclude opaque and indirect calls.}',r'\footnotesize\begin{tabular}{lrrrr}\toprule &Q ind.&T ind.&Q shared&T shared\\\midrule'])
 for label,col,data,scale in [('Host mean [us]','total_us',t,1),('Host 99th [us]','p99_us',t,1),('Jitter mean [us]','total_us',jt,1),('Target core [KiB]','sizeof_qmekf',target,1/1024),('Code+literal [KiB]','code_bytes',target,1/1024),('Read-only [KiB]','rodata_bytes',target,1/1024),('Stack-chain bound [B]','update_chain_stack',target,1)]:
  names=['A41','T41','CSM','TCSM'];rows.append(label+'&'+'&'.join(fmt(data.loc[f,col]*scale,1) for f in names)+r'\\')
 rows.extend([r'\bottomrule\end{tabular}\end{table}'])
 rows.extend([r'\begin{table*}[!t]\centering\caption{Operational consistency, averaging per-run summaries over conditions and seeds. Joint NEES is in common physical 9-state coordinates; coverage is the mean marginal 95\% position coverage.}',r'\footnotesize\begin{tabular}{lrrrrrrrr}\toprule Filter&Heave [cm]&H [cm]&V RMS [cm/s]&Tilt [deg]&Heading [deg]&NIS&NEES$_9$&Coverage [\%]\\\midrule'])
 for f in LABEL:
  g=main[main['filter']==f];vel=np.sqrt(g.rmsv_x**2+g.rmsv_y**2+g.rmsv_z**2)
  vals=[100*g.rms_z.mean(),100*g.rms_h.mean(),100*vel.mean(),g.tilt_rms_deg.mean(),g.yaw_rms_deg.mean(),g.nis.mean(),g.nees9_phys.mean(),100*g[['cov95_x','cov95_y','cov95_z']].mean().mean()]
  rows.append(tex(LABEL[f])+'&'+'&'.join(fmt(v) for v in vals)+r'\\')
 rows.extend([r'\bottomrule\end{tabular}\end{table*}'])
 rows.extend(extra_rows);rows.extend(sweep_rows)
 failures=df[df.failed>0];rows.append(r'\paragraph{Flagged trajectories.}')
 if failures.empty:rows.append('No registered estimator trajectory triggered the benchmark divergence/nonfinite flag. This flag does not test useful heading accuracy: the separately counted reduced-heave heading failures are severe practical estimation failures, even though position remains bounded.')
 else:
  for (case,f),g in failures.groupby(['case','filter']):rows.append(tex(CASE.get(case,case))+', '+tex(f)+f': {len(g)} flagged runs. ')
 psd=main.groupby('filter').agg(cov_asym_max=('cov_asym_max','max'),cov_corr_min=('cov_corr_min','min'),cov_nonfinite=('cov_nonfinite','sum'))
 psd.to_csv(results/'covariance_audit.csv')
 rows.append(r'\paragraph{Covariance audits.}')
 for f,r in psd.iterrows():rows.append(tex(LABEL[f])+f': maximum relative asymmetry {r.cov_asym_max:.3g}, minimum correlation eigenvalue {r.cov_corr_min:.3g}, nonfinite audits {int(r.cov_nonfinite)}. ')
 (DOC/'tfg-shared-tables.tex-part').write_text('\n'.join(row.rstrip() for row in rows)+'\n')
 import gzip
 with (results/'runs.csv.gz').open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
  gz.write(df.to_csv(index=False,float_format='%.9g').encode())
 manifest={'generator_sha256':__import__('hashlib').sha256(Path(__file__).read_bytes()).hexdigest(),
           'inputs':{p.name:__import__('hashlib').sha256(p.read_bytes()).hexdigest() for p in sorted(results.iterdir()) if p.is_file() and p.name not in ('analysis_manifest.json','runs.csv','reset_runs.csv') and (p.suffix in ('.csv','.json','.txt') or p.name=='runs.csv.gz')}}
 (results/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print('Article tables complete; practical improvements:',len(practical),'adverse endpoints:',len(adverse))

def main():
 """Generate publication figures and tables from committed study scores."""
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--results',type=Path,default=OUT);args=ap.parse_args();r=args.results;out=r/'figures';out.mkdir(parents=True,exist_ok=True)
 df=pd.read_csv(r/('runs.csv' if (r/'runs.csv').exists() else 'runs.csv.gz'));p=pd.read_csv(r/'paired.csv');a=pd.read_csv(r/'aggregate.csv');ab=pd.read_csv(r/'ablations.csv');inter=pd.read_csv(r/'interactions.csv')
 paired_figures(r,out,p,inter);tracking(r,out);degradation(r,out,df);ablations_and_consistency(r,out,df,ab);supplemental_figures(r,out);heave_sweep(r,out);tables_and_macros(r,df,p,a,ab)
if __name__=='__main__':main()
