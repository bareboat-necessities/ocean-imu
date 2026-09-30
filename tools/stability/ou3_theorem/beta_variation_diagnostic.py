"""Non-promoting carried diagnostic for chronological AW input weights."""
import argparse, math, subprocess, tempfile
from pathlib import Path
import numpy as np
from . import aw_tracking_source_diagnostic as base

REPO=base.REPO; DRIVER=base.DRIVER; STEP=float(base.STEP); NWIN=int(base.WINDOW_S/base.STEP)
EIGEN='/usr/include/eigen3'

ACC_OLD='''    xext.noalias() += K * r;          // State update
    joseph_update3_(K, S_mat, PCt);   // Covariance update
'''
ACC_NEW='''    if (aw_event_rec) {
        std::printf("KA %.17g %.17g %.17g %.17g", (double)qref.w(), (double)qref.x(), (double)qref.y(), (double)qref.z());
        for (int rr=0; rr<12; ++rr) for (int cc=0; cc<3; ++cc) std::printf(" %.17g", (double)K(OFF_V+rr,cc));
        std::printf("\\n");
    }
'''+ACC_OLD
S_OLD='''    xext.noalias() += K * r;            // State update
    joseph_update3_(K, S_mat, PCt);     // Covariance update
'''
S_NEW='''    if (aw_event_rec) {
        std::printf("KS");
        for (int rr=0; rr<12; ++rr) for (int cc=0; cc<3; ++cc) std::printf(" %.17g", (double)K(OFF_V+rr,cc));
        std::printf("\\n");
    }
'''+S_OLD

def compile_tapped(tmp):
    header=(REPO/base.HEADER).read_text()
    if header.count(base.TAP_ANCHOR)!=1: raise ValueError('acc row anchor changed')
    header=header.replace(base.TAP_ANCHOR,base.TAP)
    for old,new in base.EVENT_TAPS:
        if header.count(old)!=1: raise ValueError('event anchor changed')
        header=header.replace(old,new)
    if header.count(ACC_OLD)!=1 or header.count(S_OLD)!=1: raise ValueError('gain anchor changed')
    header=header.replace(ACC_OLD,ACC_NEW).replace(S_OLD,S_NEW)
    tapped=Path(tmp)/'observed/kalman_ou_iii'; tapped.mkdir(parents=True)
    (tapped/'Kalman3D_Wave_OU_III.h').write_text(header)
    exe=Path(tmp)/'observed.bin'
    subprocess.run(['g++','-O2','-std=c++20','-DAW_TRACKING_TAP','-I'+str(Path(tmp)/'observed'),
                    '-I'+str(REPO/'src'),'-I'+EIGEN,str(DRIVER),'-o',str(exe)],check=True)
    return exe

def run(exe,profile):
    name,*a=profile
    argv=[*a[:6],str(base.T0),str(base.RAMP),str(base.T_END),str(base.T0+base.RAMP),a[6]]
    lines=subprocess.run([str(exe),*argv],check=True,capture_output=True,text=True).stdout.splitlines()
    rows=[[float(x) for x in z.split()] for z in lines[:-1] if z and z[0].isdigit()]
    events=[z for z in lines[:-1] if z and not z[0].isdigit()]
    return name,rows,events

def qR(q):
    w,x,y,z=q; n=w*w+x*x+y*y+z*z
    return np.array([[w*w+x*x-y*y-z*z,2*(x*y-w*z),2*(x*z+w*y)],
      [2*(x*y+w*z),w*w-x*x+y*y-z*z,2*(y*z-w*x)],
      [2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z]])/n

def phi12(tau,h):
    e=math.exp(-h/tau); l=1/tau
    c0=(1-e)/l; c1=h/l-(1-e)/l**2; c2=h*h/(2*l)-h/l**2+(1-e)/l**3
    P=np.eye(12); I=np.eye(3)
    P[0:3,9:12]=c0*I; P[3:6,0:3]=h*I; P[3:6,9:12]=c1*I
    P[6:9,0:3]=h*h/2*I; P[6:9,3:6]=h*I; P[6:9,9:12]=c2*I; P[9:12,9:12]=e*I
    return P

def parse(events,rows):
    out=[]; ai=0; tau=rows[0][7]
    for line in events:
        f=line.split(); typ=f[0]
        if typ=='A':
            if ai<len(rows): tau=rows[ai][7]
            out.append(('A',ai)); ai+=1
        elif typ=='P': out.append(('P',phi12(tau,float(f[8]))))
        elif typ=='KA':
            K=np.array([float(x) for x in f[5:]]).reshape(12,3); R=qR([float(x) for x in f[1:5]])
            H=np.zeros((3,12)); H[:,9:12]=R
            out.append(('KA',np.eye(12)-K@H,K@R))
        elif typ=='KS':
            K=np.array([float(x) for x in f[1:]]).reshape(12,3)
            H=np.zeros((3,12)); H[:,6:9]=np.eye(3)
            out.append(('KS',np.eye(12)-K@H))
    return out

def window_bounds(ops,slide=100):
    b=np.array([float(x) for x in base.FIELD]); b/=np.linalg.norm(b)
    seed=np.array([1.,0.,0.]) if abs(b[0])<.9 else np.array([0.,1.,0.])
    u1=seed-b*np.dot(seed,b); u1/=np.linalg.norm(u1); u2=np.cross(b,u1)
    C=np.zeros((2,12)); C[0,9:12]=u1; C[1,9:12]=u2
    apos=[i for i,e in enumerate(ops) if e[0]=='A']; ans=[]
    for s in range(0,len(apos)-NWIN,slide):
        end=s+NWIN; lo,hi=apos[s],apos[end]; lam=np.zeros((12,2)); beta=[]
        for ii in range(hi,lo-1,-1):
            e=ops[ii]
            if e[0]=='A' and s<=e[1]<=end:
                lam += C.T*((.5 if e[1] in (s,end) else 1.)/NWIN)
            elif e[0]=='KA':
                A,B=e[1],e[2]; beta.append(B.T@lam); lam=A.T@lam
            elif e[0] in ('KS','P'): lam=e[1].T@lam
        beta=list(reversed(beta))
        if len(beta)<3: continue
        W=[x/STEP for x in beta]
        tv=np.linalg.norm(W[0],2)+np.linalg.norm(W[-1],2)+sum(np.linalg.norm(W[i+1]-W[i],2) for i in range(len(W)-1))
        D=[(W[i+1]-W[i])/STEP for i in range(len(W)-1)]
        tv2=np.linalg.norm(D[0],2)+np.linalg.norm(D[-1],2)+sum(np.linalg.norm(D[i+1]-D[i],2) for i in range(len(D)-1))
        ans.append((tv,tv2,float(np.linalg.norm(lam,2))))
    return ans

def main():
    global EIGEN
    ap=argparse.ArgumentParser(); ap.add_argument('--eigen',default=EIGEN); ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args(); EIGEN=a.eigen
    rec={'qualification':'OU3_BETA_VARIATION_FEASIBILITY_V1','source_uniform_certificate':False,'profiles':{}}
    with tempfile.TemporaryDirectory() as tmp:
        exe=compile_tapped(tmp)
        for p in base.PROFILES:
            if not base.admissibility(*p[1:])['admitted']: raise ValueError('inadmissible profile')
            name,rows,events=run(exe,p); vals=window_bounds(parse(events,rows))
            m1=max(x[0] for x in vals); m2=max(x[1] for x in vals); root=max(x[2] for x in vals)
            rec['profiles'][name]={'windows':len(vals),'first_variation':m1,'second_variation':m2,
              'velocity_abel_charge_mps2':float(base.LIMITS['V'])*m1,
              'displacement_raw_second_charge_mps2':float(base.LIMITS['P'])*m2,
              'root_adjoint_norm':root}
    rec['max_velocity_abel_charge_mps2']=max(x['velocity_abel_charge_mps2'] for x in rec['profiles'].values())
    rec['max_displacement_raw_second_charge_mps2']=max(x['displacement_raw_second_charge_mps2'] for x in rec['profiles'].values())
    import json; a.output.write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n'); print(json.dumps(rec,indent=2))
if __name__=='__main__': main()
