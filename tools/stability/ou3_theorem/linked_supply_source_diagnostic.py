"""80-digit linked-supply diagnostic on carried source words; no promotion.

Use --export to compile an observer and untapped control via the existing
readout fixture. Otherwise analyze previously exported word-{0,wave}.json.
Actual endpoint defects are retrospective; they are not uniform supply bounds.
"""
from pathlib import Path
import argparse
import json
import hashlib
import shutil
import subprocess
import mpmath as mp

from .ag_readout_source_diagnostic import instrument, HEADER
REPO = Path(__file__).resolve().parents[3]

def mat(a): return mp.matrix([[mp.mpf(float(x)) for x in row] for row in a])
def sym(a): return (a+a.T)/2
def skew(v):
 x,y,z=v
 return mp.matrix([[0,-z,y],[z,0,-x],[-y,x,0]])
def qprod(a,b):
 va,vb=mp.matrix(a[:3]),mp.matrix(b[:3]); v=a[3]*vb+b[3]*va+skew(va)*vb
 return list(v)+[a[3]*b[3]-(va.T*vb)[0]]
def error(state, qraw, t, live, wave):
 # NED, error q_true=Exp(delta_theta)q_hat (left W->B error).
 q=[mp.mpf(float(row[0])) for row in qraw]; norm=mp.sqrt(sum(x*x for x in q));q=[x/norm for x in q]
 angle=-mp.mpf('.02')*mp.sin(t/2) if wave else mp.mpf(0)
 qt=[mp.sin(angle/2),mp.mpf(0),mp.mpf(0),mp.cos(angle/2)]
 qe=qprod(qt,[-q[0],-q[1],-q[2],q[3]])
 if qe[3]<0: qe=[-x for x in qe]
 n=mp.sqrt(sum(x*x for x in qe[:3]));rv=mp.matrix(qe[:3])*(2*mp.atan2(n,qe[3])/n) if n else mp.zeros(3,1)
 e=-mat(state);e[:3,0]=rv
 if wave:
  e[8]+=mp.mpf('.24')*mp.cos(mp.mpf('.6')*t)
  e[11]+=mp.mpf('.4')*mp.sin(mp.mpf('.6')*t)
  e[14]+=mp.mpf(2)/3*(mp.cos(mp.mpf('.6')*live)-mp.cos(mp.mpf('.6')*t))
  e[17]-=mp.mpf('.144')*mp.sin(mp.mpf('.6')*t)
 return e
fmt=lambda x:mp.nstr(x,25)
def analyze(path):
 trace=json.loads(path.read_text()); wave=trace['profile']=='wave'
 M=mp.eye(21)
 for ev in trace['events']:
  if ev['kind']=='prediction':
   M[:6,:]=mat(ev['F_AG'])*M[:6,:]
   M[6:18,:]=mat(ev['F_LIN'])*M[6:18,:]
   M[18:,:]=mp.mpf(float(ev['phi_BA']))*M[18:,:]
  elif ev['kind']=='correction': M=M-mat(ev['K'])*(mat(ev['H'])*M)
  elif ev['kind']=='reset':
   v=[mp.mpf(float(x[0])) for x in ev['d']]
   M[:3,:]=(mp.eye(3)+skew(v)/2)*M[:3,:]
  elif ev['kind'] not in ('sync','sync_completion'): raise ValueError(ev['kind'])
 P0=sym(mat(trace['root_covariance'])); PN=sym(mat(trace['terminal_covariance']))
 L0=mp.cholesky(P0);LN=mp.cholesky(PN); Li=LN**-1
 H=Li*M*L0
 Delta=sym(mp.eye(21)-H.T*H); eig=mp.eigsy(Delta,eigvals_only=True); dmin=min(eig)
 t0=mp.mpf(225);tn=mp.mpf('225.32');live=mp.mpf(trace['live_step'])/200
 e0=error(trace['root_state'],trace['root_quaternion'],t0,live,wave)
 e1=error(trace['terminal_state'],trace['terminal_quaternion'],tn,live,wave)
 x=L0**-1*e0; b=Li*(e1-M*e0); z=H.T*b
 V0=(x.T*x)[0]; VN=(e1.T*(PN**-1)*e1)[0]
 ans={'profile':trace['profile'],'dps':80,'event_count':len(trace['events']),
  'observer_terminal_parity':trace['observer_terminal_parity'],
  'trace_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
  'source_header_sha256':trace['source_header_sha256'],
  'root_sqrt_V':fmt(mp.sqrt(V0)),'terminal_sqrt_V':fmt(mp.sqrt(VN)),
  'V_change':fmt(VN-V0),'homogeneous_rho':fmt(1-dmin),
  'minimum_whitened_loss':fmt(dmin),'signed_endpoint_forcing_norm':fmt(mp.norm(b)),
  'active_release_time_s':trace['active_step']/200,
  'live_time_s':trace['live_step']/200,
  'release_LIN_mean_norm':fmt(mp.norm(mat(trace['release_state'])[6:18,:])),
  'source_uniform_verified':False,'theorem_closed':False,
  'supply_is_retrospective_endpoint_residual':True,
  'uniform_inner_retention_refuted_by_this_finite_replay':False,
  'all_time_magnetic_service_verified':False,
  'physical_S_origin':'one fixed capture/Live epoch; never reset at word boundaries',
  'metric':'exact symmetric part of exported float covariance; arithmetic enclosure not claimed'}
 if dmin<=0: ans['status']='NO_POSITIVE_FROZEN_WORD_LOSS';return ans
 def eval_gamma(g):
  GG=Delta-g*mp.eye(21)
  cc=(b.T*b)[0]+(z.T*(GG**-1)*z)[0]
  return cc/g,cc,GG
 # Convex-looking one-dimensional objective; use a dense logit grid then
 # golden refinement inside the best bracket. Diagnostic only.
 grid=[dmin*mp.mpf(k)/1000 for k in range(1,1000)]
 vals=[eval_gamma(g)[0] for g in grid]
 ib=min(range(len(vals)),key=lambda k: vals[k])
 lo=grid[max(0,ib-1)]; hi=grid[min(len(grid)-1,ib+1)]
 gr=(mp.sqrt(5)-1)/2
 a,c=lo,hi; x1=c-gr*(c-a); x2=a+gr*(c-a); f1=eval_gamma(x1)[0]; f2=eval_gamma(x2)[0]
 for _ in range(120):
  if f1>f2: a=x1;x1=x2;f1=f2;x2=a+gr*(c-a);f2=eval_gamma(x2)[0]
  else: c=x2;x2=x1;f2=f1;x1=c-gr*(c-a);f1=eval_gamma(x1)[0]
 gam=(a+c)/2; ratio,chi,G=eval_gamma(gam)
 completed=(1-gam)*V0+chi-((x-G**-1*z).T*G*(x-G**-1*z))[0]
 ans.update({'defect_origin':'retrospective_endpoint_residual',
  'eligible_for_theorem_entry_test':False,
  'gamma':fmt(gam),'linked_supply_chi':fmt(chi),
  'optimized_chi_over_gamma':fmt(ratio),
  'optimized_linked_to_inner_budget_ratio':fmt(ratio/(mp.mpf('.15')**2)),
  'frozen_fixed_forcing_sufficient_radius':fmt(mp.sqrt(ratio)),
  'linked_to_inner_budget_ratio':fmt(chi/(gam*mp.mpf('.15')**2)),
  'linked_to_root_budget_ratio':fmt(chi/(gam*V0)) if V0 else None,
  'separated_fixed_forcing_radius':fmt(mp.norm(b)/(1-mp.sqrt(1-dmin))),
  'completed_square_identity_error':fmt(abs(VN-completed)),
  'status':'NON_PROMOTING_LINKED_FEASIBILITY_ONLY'})
 return ans


def export_words(directory, eigen, compiler):
    directory.mkdir(parents=True, exist_ok=True)
    inc = directory/'include'/'kalman_ou_iii'
    inc.mkdir(parents=True, exist_ok=True)
    header = (REPO/HEADER).read_text()
    (inc/HEADER.name).write_text(instrument(header))
    source = (REPO/'tools/stability/ag_readout_source.cpp').read_text()
    def once(old, new):
        nonlocal source
        if source.count(old) != 1:
            raise ValueError('readout driver anchor changed: '+old[:70])
        source = source.replace(old, new)
    once('    std::string root;', '    std::string root, root_state, root_quat, release_state;')
    once('            recording=true;',
         '            root_state=matrix_json(filter.raw().mekf().xext);\n'
         '            root_quat=matrix_json(filter.raw().mekf().qref.coeffs());\n'
         '            recording=true;')
    once('if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;',
         'if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) {\n'
         '          active=k; release_state=matrix_json(filter.raw().mekf().xext);\n        }')
    needle = '              << ",\\\"terminal_covariance\\\":" << matrix_json(m.Pext)'
    # Avoid altering any estimator or shared driver in the repository.
    once(needle,
         '              << ",\\\"root_state\\\":" << root_state\n'
         '              << ",\\\"root_quaternion\\\":" << root_quat\n'
         '              << ",\\\"release_state\\\":" << release_state\n'+needle)
    driver = directory/'linked-source.cpp'
    driver.write_text(source)
    source_hashes = {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted((REPO/'src').rglob('*.h'))}
    for name, includes in [('observed', ['-I'+str(inc.parent)]), ('control', [])]:
        subprocess.run([compiler, '-O1', '-std=c++20', *includes,
                        '-I'+str(REPO/'src'), '-isystem', str(eigen),
                        str(driver), '-o', str(directory/name)], check=True)
    for profile in ['0', 'wave']:
        obs = json.loads(subprocess.check_output([str(directory/'observed'), profile], text=True))
        ctl = json.loads(subprocess.check_output([str(directory/'control'), profile], text=True))
        if any(obs[k] != v for k, v in ctl.items() if k != 'events'):
            raise ArithmeticError('observer changed source execution')
        obs.update(profile=profile, observer_terminal_parity=True,
                   source_header_sha256=hashlib.sha256(header.encode()).hexdigest(),
                   source_hashes=source_hashes)
        (directory/('word-'+profile+'.json')).write_text(json.dumps(obs, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--export', action='store_true')
    parser.add_argument('--eigen', type=Path, default=Path('/usr/include/eigen3'))
    parser.add_argument('--cxx', default=shutil.which('clang++') or 'g++')
    args = parser.parse_args()
    if args.export:
        export_words(args.directory.resolve(), args.eigen, args.cxx)
    with mp.workdps(80):
        report = {'kind': 'NON_PROMOTING_LINKED_FINITE_SUPPLY', 'decimal_digits': 80,
                  'source_uniform_verified': False, 'theorem_closed': False,
                  'words': [analyze(args.directory/('word-'+p+'.json')) for p in ['0', 'wave']]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
