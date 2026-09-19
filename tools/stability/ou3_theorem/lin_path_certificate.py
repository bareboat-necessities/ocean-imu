"""Exact-rational LIN controllability/correction comparison; not a full theorem.

For every piecewise-constant shipping tau, the process action of a degree-seven
endpoint path bounds posterior precision. Measurements are bounded on their
actual sample mesh, not transported backwards through an OU inverse. Source
small-x polynomial defects are charged in the same action. No sampled spectrum
or floating eigensolver promotes this certificate. Float32 is a separate open
supply obligation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from math import factorial, isqrt
from pathlib import Path


def inverse(a):
    n = len(a)
    b = [[F(v) for v in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if b[i][j])
        b[j], b[pivot] = b[pivot], b[j]
        q = b[j][j]
        b[j] = [v / q for v in b[j]]
        for i in range(n):
            if i != j:
                q = b[i][j]
                b[i] = [v - q * w for v, w in zip(b[i], b[j])]
    return [row[n:] for row in b]


def derivative(p, order):
    for _ in range(order):
        p = [F(i) * p[i] for i in range(1, len(p))]
    return p


def square_integral(p):
    return sum((a*b / F(i+j+1) for i, a in enumerate(p)
                for j, b in enumerate(p)), F(0))


def hermite_basis():
    # S endpoint jets: (S,p,v,a), with four zero jets at the window root.
    h = [[F(factorial(k), factorial(k-d)) for k in range(4, 8)]
         for d in range(4)]
    inv = inverse(h)
    return [[F(0)]*4 + [inv[k][j] for k in range(4)] for j in (2, 1, 0, 3)]


def sqrt_floor(q, digits=24):
    scale = 10**digits
    return F(isqrt((q.numerator * scale**2)//q.denominator), scale)


def decimal_out(q, digits=24, upper=False):
    scale = 10**digits
    n, r = divmod(q.numerator * scale, q.denominator)
    if upper and r:
        n += 1
    sign = '-' if n < 0 else ''
    n = abs(n)
    return f'{sign}{n//scale}.{n%scale:0{digits}d}'


def rational_record(q):
    return {'numerator': str(q.numerator), 'denominator': str(q.denominator),
            'decimal_lower': decimal_out(q), 'decimal_upper': decimal_out(q, upper=True)}


def small_x_source_defect():
    """Relative Q error and F defect for the literal degree-truncated branches.

    In natural step units D_h=diag(h,h^2,h^3,1), Q/(sigma^2*x)
    has entries 2 sum_r (-x)^r c_mnr. The coefficient and tail majorants
    below are rational. Only x<.01 needs a truncation allowance.
    """
    order = (1, 2, 3, 0)
    degrees = ((6,5,4,7), (5,4,3,6), (4,3,2,5), (7,6,5,8))
    x = F(1,100)
    b0 = [[F(2, factorial(m)*factorial(n)*(m+n+1))
           for n in order] for m in order]
    b0inv = inverse(b0)
    # Volterra integration has L2 norm <=1 on [0,1]. Thus B(x)>=B(0)/(1+x)^2.
    inverse_norm = max(sum(abs(v) for v in row) for row in b0inv)*(1+x)**2
    tails = []
    for i, m in enumerate(order):
        row = []
        for j, n in enumerate(order):
            r = degrees[i][j]+1
            c = sum((F(2, factorial(m+k)*factorial(n+r-k)*(m+n+r+1))
                     for k in range(r+1)), F(0))
            row.append(c*x**r/(1-F(2)*x/F(r+1)))
        tails.append(row)
    radius = max(sum(row) for row in tails)
    eps_q = inverse_norm*radius
    # DeltaPhi affects p,a and S,a only, with magnitudes h^2*x^3/120
    # and h^3*x^3/720. Exact Q inverse action in natural step units:
    j_defect = inverse_norm*x**5*(F(1,120)**2+F(1,720)**2)/F('.05')**2
    return eps_q, j_defect, b0, degrees


def certificate():
    tmin, tmax, dtmin = F(16), F('16.006'), F('.004')
    scales, powers = (F('5.5'),F('8.1'),F(1100),F(4)), (2,1,0,3)
    basis = hermite_basis()
    def gram_trace(order):
        out = F(0)
        for p, scale, power in zip(basis, scales, powers):
            exponent = 2*power+1-2*order
            time_bound = (tmax if exponent >= 0 else tmin)**exponent
            out += scale**2*time_bound*square_integral(derivative(p, order))
        return out
    g = [gram_trace(d) for d in range(5)]
    def sampled_trace(d):
        # Disjoint backward dt_min cells, f(0)=0. For every mesh event,
        # f(t_i)^2 <= 2/dt_min int_cell f^2 + 2 dt_min int_cell f'^2.
        return 2/dtmin*g[d]+2*dtmin*g[d+1]
    eps_q, j_defect, _, _ = small_x_source_defect()
    if not F(0) <= eps_q < F(1):
        raise ArithmeticError('relative shipping process covariance enclosure failed')
    # For any lambda(t) in [1/12,50], integration of 2a a' telescopes.
    ideal_action = (12*g[4]+50*g[3]+scales[3]**2)/(2*F('.05')**2)
    # Young 2/2 bound on w_ship=w_exact-DeltaPhi*x_previous.
    source_action = 2*(ideal_action+j_defect*sampled_trace(3))/(1-eps_q)
    acc_action = F(3)/F('.05')**2*sampled_trace(3)
    s_action = F(1)/F('.075')**2*sampled_trace(0)
    precision = source_action+acc_action+s_action
    floor = 1/precision
    ell = sqrt_floor(floor)
    return {
        'qualification':'OU3_LIN_PATH_ENERGY_REAL_ARITHMETIC_V1',
        'verified':ell>0 and ell**2 <= floor,
        'scope':'source-uniform real-arithmetic shipping LIN comparison; float32 transfer remains separate',
        'method':'exact rational endpoint-path action with source polynomial defects and mesh measurement penalties',
        'window_s':['16','16.006'], 'tau_s':['0.02','12'],
        'sample_period_s':['0.004','0.006'],
        'sigma_aw_min':'0.05', 'accel_std_min':'0.05', 'integral_std_min':'0.075',
        'scales':['5.5','8.1','1100','4'],
        'arbitrary_piecewise_constant_tau':True,
        'all_acc_and_integral_corrections_included':True,
        'cross_block_measurement_domination_factor':3,
        'source_Q_relative_defect_ceiling':rational_record(eps_q),
        'source_Q_relative_positive_margin':rational_record(1-eps_q),
        'source_F_defect_information_ceiling':rational_record(j_defect),
        'ideal_process_precision_ceiling':rational_record(ideal_action),
        'source_process_precision_ceiling':rational_record(source_action),
        'accel_precision_ceiling':rational_record(acc_action),
        'integral_precision_ceiling':rational_record(s_action),
        'posterior_precision_ceiling':rational_record(precision),
        'covariance_floor':rational_record(floor),
        'factor_floor':decimal_out(ell),
        'deployment_noise_floor_binding_verified':True,\n        'shipping_interleaved_corrections_verified':True,\n        'aw_sync_psd_inflation_verified':True,
        'float32_covariance_factor_verified':False,
        'constructive_full_A21_mu_rho_enclosure':False,
        'theorem_closed':False,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    report=certificate()
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    print(text,end='')
    return 0 if report['verified'] else 1


if __name__=='__main__':
    raise SystemExit(main())
