#!/usr/bin/env python3
"""Generate the LaTeX evidence tables and figures from the verified replay rows."""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t as student

from marine_mag_study import OUT,ROOT,NAV,plan
sys.path.insert(0, str(ROOT / 'plots'))
from vector_export import save_figure

LABEL={-1:'Sketch opt-out',0:'Replay default',1:'Existing enabled',2:'Hard iron',3:'Diagonal',
       4:'Symmetric',5:'Full affine',6:'Selective',7:'Affine + covariance',8:'V7 + quarantine',9:'V6 + covariance'}
DOC=ROOT/'doc/marine_mag'
MATCH=['case','spectrum','height','seed','scenario','dip']


def paired(d,metric):
    base=d[d.variant==0].set_index(MATCH)[metric].rename('baseline')
    return d.join(base,on=MATCH).assign(delta=lambda x:x[metric]-x.baseline)


def statistics(d):
    rows=[]
    for metric in [*NAV,'centered_deg','yaw_mean_deg','matrix_err_rms','applied_matrix_err_rms']:
        for v,g in paired(d,metric).groupby('variant'):
            draws=g.groupby('seed').delta.mean()
            h=student.ppf(.975,len(draws)-1)*draws.std()/np.sqrt(len(draws))
            rows.append(dict(metric=metric,variant=v,mean=g[metric].mean(),delta=g.delta.mean(),
                             pct=100*g.delta.mean()/g.baseline.mean(),ci_low=draws.mean()-h,ci_high=draws.mean()+h,
                             improved=sum(g.delta < -1e-6),tied=sum(g.delta.abs()<=1e-6),worse=sum(g.delta>1e-6),
                             p95=g[metric].quantile(.95),worst=g[metric].max(),worst_delta=g.delta.max()))
    return pd.DataFrame(rows)


def table(file,rows):
    file.write_text('\n'.join(' & '.join(map(str,r))+r' \\' for r in rows)+'\n'+r'\bottomrule'+'\n')


def figures(info,traces,stats,doc):
    # Reuse the repository's figure envelope/decimation and embedded-font settings.
    plt.rcParams.update({'font.family':'serif','font.size':9,'svg.fonttype':'none','pdf.fonttype':42})
    fig,ax=plt.subplots(1,2,figsize=(7,2.35),sharey=True,layout='constrained')
    for a,case in zip(ax,('following','head')):
        s=info[(info.case==case)&(info.height==1.5)&(info.window==1200)&(info.waves==1)]
        for model,label in [('hard','Offsets'),('diagonal','Offsets + diagonal'),('affine','Scale-fixed affine')]:
            d=s[(s.model==model)&(s.steering!='turn')].copy()
            d['amp']=d.steering.astype(int)
            d=d.sort_values('amp')
            a.semilogy(d.amp,d.smin,'o-',lw=1,ms=3,label=label)
        a.axhline(2,color='.4',ls=':',lw=1)
        a.set(title=case.capitalize()+' seas',xlabel='Steering amplitude (degrees)',xticks=[0,1,3,5,10])
    ax[0].set_ylabel('Weakest singular value')
    ax[1].legend(fontsize=7)
    save_figure(fig,doc/'observability.pdf',metadata={'CreationDate':None});plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(6.9,3.5),sharex=True,layout='constrained')
    for v in [0,6,7]:
        d=traces[traces.variant==v]
        axes[0].plot(d.t,d.yaw_mean,label=LABEL[v],lw=1)
        axes[1].plot(d.t,np.linalg.norm(d[['bx','by','bz']],axis=1),label=LABEL[v],lw=1)
    axes[0].set_ylabel('Signed yaw error (deg)')
    axes[1].set_ylabel('Offset error norm (µT)')
    axes[1].set_xlabel('Elapsed time (s)')
    axes[0].legend(ncol=3,fontsize=8)
    save_figure(fig,doc/'convergence.pdf',metadata={'CreationDate':None});plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.8,2.7),layout='constrained')
    d=stats[(stats.metric=='yaw_rms_deg')&stats.variant.isin([-1,2,3,4,6,7,8])].sort_values('variant')
    ax.errorbar(d.delta,np.arange(len(d)),xerr=[d.delta-d.ci_low,d.ci_high-d.delta],
                fmt='o',color='#145f73',capsize=3)
    ax.axvline(0,ls=':',lw=.8,color='.3');ax.set_yticks(np.arange(len(d)),[LABEL[v] for v in d.variant])
    ax.set_xlabel('Paired change in yaw RMS vs V0 (degrees), 95% interval')
    save_figure(fig,doc/'paired.pdf',metadata={'CreationDate':None});plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=OUT)
    p.add_argument('--doc',type=Path,default=DOC)
    args=p.parse_args();o=args.out;doc=args.doc
    doc.mkdir(parents=True,exist_ok=True)
    for stage in ['diagnosis','development','validation','stress']:
        frame=pd.read_csv(o/(stage+'.csv'))
        assert len(frame)==len(list(plan(stage))) and not frame.duplicated(MATCH+['variant']).any(), stage
        # V0/V1/V5 are deliberately separate executions of the same model.
        if stage!='diagnosis':
            for m in NAV:
                piv=frame.pivot(index=MATCH,columns='variant',values=m)
                assert np.allclose(piv[0],piv[1],atol=1e-7,rtol=0) and np.allclose(piv[0],piv[5],atol=1e-7,rtol=0),m
    d=pd.read_csv(o/'validation.csv')
    st=statistics(d); st.to_csv(o/'paired.csv',index=False,float_format='%.8g')
    quarantine=pd.read_csv(o/'quarantine-only.csv').assign(variant=10)
    trial=pd.concat([d[d.variant==0],quarantine],ignore_index=True)
    ablation=statistics(trial)
    ablation[ablation.variant==10].to_csv(o/'quarantine-paired.csv',index=False,float_format='%.8g')
    y=st[st.metric=='yaw_rms_deg'].set_index('variant')
    table(doc/'paired-rows.tex',[(LABEL[v],f'{r["mean"]:.3f}',f'{r.delta:+.3f}',f'{r.pct:+.1f}',
          f'[{r.ci_low:+.2f}, {r.ci_high:+.2f}]',f'{int(r.improved)}/{int(r.worse)}',f'{r.worst_delta:+.2f}')
          for v,r in y.iterrows()])
    nav=d.groupby('variant')[['yaw_mean_deg','centered_deg','roll_rms_deg','pitch_rms_deg','disp_3d_rms_m','gyro_bias_3d_rms_radps','accel_bias_3d_rms_mps2']].mean()
    table(doc/'nav-rows.tex',[(v,*[f'{n:.3f}' if i!=5 else f'{1e3*n:.3f}' for i,n in enumerate(row)]) for v,row in nav.iterrows()])
    q=pd.read_csv(o/'diagnosis.csv')
    table(doc/'cause-rows.tex',[(s.replace('_',r'\_'),*[f'{q[(q.scenario==s)&(q.case==case)&(q.variant==v)].yaw_rms_deg.iloc[0]:.3f}'
           for case in ('following','head') for v in (-1,0)]) for s in q.scenario.unique()])
    figures(pd.read_csv(o/'information.csv'),pd.read_csv(o/'development-trace.csv'),st,doc)


if __name__=='__main__':
    main()
