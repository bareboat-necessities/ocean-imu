"""Optional 80-digit supplied-word feasibility probe; never promotes a claim.

Run python -m tools.stability.ou3_theorem.ag_readout_diagnostic (needs mpmath).
Exact enclosure of the conditional action uses ag_readout.certificate;
neither this test word nor its prior family is claimed shipping-reachable.
"""
import json

from .ag_readout import exact_readout, readout_action, supplied_fixture
from .matrix_certificates import identity


def diagnostic(dps=80):
    import mpmath as mp

    if dps < 80:
        raise ValueError('use at least 80 decimal digits')
    with mp.workdps(dps):
        def matrix(a):
            return mp.matrix([[mp.mpf(v.numerator)/v.denominator for v in row] for row in a])

        events = supplied_fixture()
        reader = exact_readout(events)
        audit = readout_action(events, reader, identity(15))
        bound = matrix(audit['action'])
        results = []
        for scale in (1, 10**6, 10**12):
            c = mp.eye(21)
            for i in range(6):
                c[i, i] = mp.sqrt(scale)
                c[i, 6+i] = mp.mpf(1)/4
            p, root_factor = c*c.T, c.copy()
            root = p.copy()
            transport = mp.eye(21)
            for event in events:
                kind = event['kind']
                if kind in ('prediction', 'reset'):
                    f = matrix(event['F'] if kind == 'prediction' else event['G'])
                    p = f*p*f.T
                    transport = f*transport
                    if kind == 'prediction':
                        u = matrix(event['U'])
                        p += u*u.T
                else:
                    h, v = matrix(event['H']), matrix(event['V'])
                    r = v*v.T
                    k = p*h.T*(h*p*h.T+r)**-1
                    a = mp.eye(21)-k*h
                    p = a*p*a.T+k*r*k.T
                    transport = a*transport
            endpoint = transport.T*(p**-1)*transport
            loss = root**-1-endpoint  # high precision diagnostic only
            rho = max(mp.eigsy(root_factor.T*endpoint*root_factor, eigvals_only=True))
            results.append({
                'root_AG_scale': str(scale),
                'AG_trial_slack_min': mp.nstr(min(mp.eigsy(bound-p[:6, :6], eigvals_only=True)), 24),
                'AG_loss_min': mp.nstr(min(mp.eigsy(loss[:6, :6], eigvals_only=True)), 24),
                'rho': mp.nstr(rho, 24)})
        return {'qualification': 'OU3_AG_READOUT_FEASIBILITY_ONLY_V1', 'decimal_digits': dps,
                'fixture': 'varying rational 21-state supplied word, correlated process factors and nonorthogonal resets',
                'AG_trial_ceiling_max': mp.nstr(max(mp.eigsy(bound, eigvals_only=True)), 24),
                'cases': results, 'source_uniform_verified': False,
                'shipping_history': False, 'theorem_closed': False}


if __name__ == '__main__':
    print(json.dumps(diagnostic(), indent=2, sort_keys=True))
