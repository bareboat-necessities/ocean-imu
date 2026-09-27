"""Carried-word adjoint feasibility, not a source-uniform separation certificate.

Reuse the historical-reader observer on temporary copies of shipping source.
Check terminal parity against its untapped control. No reseed or setter is
introduced. The four actual S events define the existing signed spline.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from .ag_readout_source_diagnostic import HEADER, REPO, instrument
from .matrix_certificates import identity, matmul, encoded
from .signed_temporal import adjoint_compatibility, atoms, jet


def driver_source():
    source=(REPO/'tools/stability/ag_readout_source.cpp').read_text()
    # One second supplies four applied S events under the literal scheduler.
    for old,new in [('k<=45064','k<=45200'),('applied!=8','applied!=25')]:
        if source.count(old)!=1: raise ValueError('readout driver anchor changed')
        source=source.replace(old,new)
    return source


def carried_functional(trace):
    """Return exact-rational C,W for the full 18 Euclidean mean coordinates.

    C_i=Phi_(N,i+1) K_i, with prediction LIN/BA transports, BG identity,
    additive correction jumps and no Euclidean change at attitude reset.
    Quaternion dependence is in actual K and W; resets/arithmetic/projection
    remain in the literal affine defect, never in an invented sensor noise.
    """
    timed=[]; step=0
    for event in trace['events']:
        if event['kind']=='prediction': step+=1
        timed.append((F(step,200),event))
    knots=[t for t,e in timed if e['kind']=='correction' and e['sensor']=='S'][:4]
    if len(knots)!=4: raise ValueError('four actual S events required')
    word=[(t,e) for t,e in timed if t<=knots[-1]]
    s_sites={t:i for i,(t,e) in enumerate(word)
             if e['kind']=='correction' and e['sensor']=='S' and t in knots}
    c_atom=dict(zip(knots,atoms(knots)))
    transport=identity(18); columns=[]; targets=[]; labels=[]
    for index in reversed(range(len(word))):
        t,e=word[index]
        if e['kind']=='prediction':
            a=identity(18)
            for i in range(12):
                for j in range(12): a[3+i][3+j]=F(e['F_LIN'][i][j])
            for i in range(15,18): a[i][i]=F(e['phi_BA'])
            transport=matmul(transport,a)
        elif e['kind']=='correction':
            # psi'' jumps at the S atom. Other same-time corrections before
            # that atom use its LEFT value; S itself and later events use RIGHT.
            side='left' if t in s_sites and index<s_sites[t] else 'right'
            k=[[F(x) for x in row] for row in e['K']]
            c=matmul(transport,k[3:])
            w=[[-jet(knots,t)*k[6+i][j]+jet(knots,t,1)*k[9+i][j]
                -jet(knots,t,2,side)*k[12+i][j]
                +(c_atom.get(t,F(0)) if e['sensor']=='S' and i==j else 0)
                for j in range(3)] for i in range(3)]
            columns[:0]=[list(col) for col in zip(*c)]
            targets[:0]=[list(col) for col in zip(*w)]
            labels[:0]=[{'t':str(t),'sensor':e['sensor'],'axis':j} for j in range(3)]
    return [list(r) for r in zip(*columns)], [list(r) for r in zip(*targets)], labels, knots


def analyze(trace):
    import mpmath as mp
    c,w,labels,knots=carried_functional(trace)
    # Quiet source has invariant coordinate channels. Restrict columns of the
    # x functional (acc_x, S_x, mag_z) WITHOUT deleting any nonzero state row.
    # Any incompatibility on this subset disproves compatibility on the word.
    indices=[j for j,l in enumerate(labels) if
             (l['sensor'] in ('acc','S') and l['axis']==0) or
             (l['sensor']=='mag' and l['axis']==2)]
    rows=[i for i in range(18) if any(c[i][j] for j in indices)]
    cs=[[c[i][j] for j in indices] for i in rows]
    ws=[[w[0][j] for j in indices]]
    # Exact elimination supplies a kernel witness; high precision measures
    # its size. No subdivision/enclosure of a failed compatibility is attempted.
    result=adjoint_compatibility(cs,ws)
    report={'scope':'one exported source word interpreted as exact rational coefficients',
            'knots_s':[str(x) for x in knots], 'state_rows':[i+3 for i in rows],
            'source_uniform_verified':False,'theorem_closed':False,
            'compatible':result['compatible']}
    if not result['compatible']:
        v=result['kernel_witness']; scale=max(abs(x) for x in v)
        v=[x/scale for x in v]
        used=[i for i,x in enumerate(v) if x]
        cf=[[r[i] for i in used] for r in cs]; wf=[[r[i] for i in used] for r in ws]
        vf=[[v[i]] for i in used]; residual=matmul(wf,vf)[0][0]
        assert matmul(cf,vf)==[[0] for _ in cf]
        assert residual!=0
        with mp.workdps(80):
            number=lambda x: mp.mpf(x.numerator)/x.denominator
            value=abs(number(residual))
            lower=F(str(mp.nstr(value/2,12)))
            assert 0<lower<abs(residual)
            report.update({'normalized_residual_80_digits':mp.nstr(value,32),
                           'residual_strict_lower':str(lower),
                           'C_kernel_exact':True,'W_kernel_nonzero_exact':True,
                           'selected_columns':[labels[indices[i]] for i in used],
                           'C':encoded(cf),'W':encoded(wf),'kernel_witness':encoded(vf)})
    return report


def run(eigen):
    source=(REPO/HEADER).read_text(); driver=driver_source()
    with tempfile.TemporaryDirectory(prefix='ou3-signed-') as directory:
        tmp=Path(directory); include=tmp/'kalman_ou_iii'; include.mkdir()
        (include/HEADER.name).write_text(instrument(source))
        cpp=tmp/'observer.cpp'; cpp.write_text(driver)
        records=[]
        for name,inc in [('observed',['-I'+str(tmp)]),('control',[])]:
            binary=tmp/name
            subprocess.run(['g++','-O2','-std=c++20',*inc,'-I'+str(REPO/'src'),
                            '-isystem',str(eigen),str(cpp),'-o',str(binary)],check=True)
            records.append(json.loads(subprocess.check_output([str(binary),'0'],text=True)))
        observed,control=records
        for key in control:
            if key!='events' and observed[key]!=control[key]:
                raise ValueError('observer changed literal terminal output: '+key)
        return {'qualification':'OU3_CARRIED_SIGNED_ADJOINT_DIAGNOSTIC_V1',
                'shipping_header_sha256':hashlib.sha256(source.encode()).hexdigest(),
                'instrumented_header_sha256':hashlib.sha256(instrument(source).encode()).hexdigest(),
                'driver_sha256':hashlib.sha256(driver.encode()).hexdigest(),
                'trace_sha256':hashlib.sha256(json.dumps(observed,sort_keys=True).encode()).hexdigest(),
                'literal_terminal_parity':True,'decimal_digits':80,
                'profile':'quiet truth; construction through release; 225 to 226 seconds',
                'all_time_magnetic_service_certified':False,
                **analyze(observed)}


def verify_diagnostic(report):
    """Independently verify the committed finite-word incompatibility witness."""
    c=[[F(x) for x in row] for row in report['C']]
    w=[[F(x) for x in row] for row in report['W']]
    v=[[F(x) for x in row] for row in report['kernel_witness']]
    if (report['compatible'] or not report['literal_terminal_parity']
            or report['source_uniform_verified'] or report['theorem_closed']):
        raise ValueError('invalid diagnostic scope')
    if max(abs(row[0]) for row in v)!=1:
        raise ValueError('kernel witness normalization changed')
    if matmul(c,v)!=[[0] for _ in c]: raise ValueError('C kernel cancellation failed')
    lower=F(report['residual_strict_lower'])
    if not 0<lower<abs(matmul(w,v)[0][0]):
        raise ValueError('W kernel separation failed')
    if report['driver_sha256']!=hashlib.sha256(driver_source().encode()).hexdigest():
        raise ValueError('signed observer driver fingerprint changed')
    if report['shipping_header_sha256']!=hashlib.sha256((REPO/HEADER).read_bytes()).hexdigest():
        raise ValueError('signed observer shipping fingerprint changed')
    if report['instrumented_header_sha256']!=hashlib.sha256(instrument((REPO/HEADER).read_text()).encode()).hexdigest():
        raise ValueError('signed observer instrumentation fingerprint changed')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen',type=Path,default=Path('/usr/include/eigen3'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(run(args.eigen),indent=2,sort_keys=True)+'\n')
