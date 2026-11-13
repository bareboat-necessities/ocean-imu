#!/usr/bin/env python3
"""Independent exact colored-noise covariance, with measurement feedback."""
import numpy as np
for variable,fixed_exposure in ((False,False),(True,False),(True,True)):
    joint=np.diag([0.,1.]);P=naive=0.;last_dt=.1;worst=0
    for k in range(1000):
        dt=.1*(3 if variable and k%17==0 else .7 if variable and k%2 else 1.)
        a=.95;gain=1-a;Q=.1;R=Q/(.1 if fixed_exposure else dt);rv=.25
        A=np.array([[a,a*.5*dt],[0,0]]);B=np.array([a*.5*dt,1]);D=np.array([gain,0])
        joint=A@joint@A.T+R*np.outer(B,B)+rv*np.outer(D,D)
        scale=max(1.,dt/.1) if fixed_exposure else max(1.,.5*(1+dt/last_dt));last_dt=dt
        P=a*a*(P+Q*dt*scale)+gain*gain*rv
        naive=a*a*(naive+.5*Q*dt)+gain*gain*rv
        worst=max(worst,joint[0,0]/P)
        assert P+1e-12>=joint[0,0]
    assert joint[0,0]>naive*1.3
    print(dict(variable_dt=variable,fixed_exposure=fixed_exposure,largest_actual_to_bound=worst,actual_to_naive_half_noise=joint[0,0]/naive))
