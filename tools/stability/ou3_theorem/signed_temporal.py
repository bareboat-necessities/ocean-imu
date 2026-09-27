"""Exact fail-closed checks for the OU-III signed temporal construction."""
from __future__ import annotations
from fractions import Fraction as F
from math import prod

def atoms(knots):
    k=tuple(F(x) for x in knots)
    if len(k)!=4 or any(b<=a for a,b in zip(k,k[1:])):
        raise ValueError("four increasing actual-S times required")
    return tuple(-F(6,1)/prod(t-s for i,s in enumerate(k) if i!=j) for j,t in enumerate(k))

def jet(knots,t,order=0,side="right"):
    k=tuple(F(x) for x in knots); t=F(t); c=atoms(k)
    if order not in (0,1,2): raise ValueError("order")
    ans=F(0)
    for cj,tj in zip(c,k):
        active=t>tj or (t==tj and side=="right")
        if active:
            if order==0: ans += cj*(t-tj)**2/F(2)
            elif order==1: ans += cj*(t-tj)
            else: ans += cj
    return ans

def endpoint_jets_vanish(knots):
    k=tuple(F(x) for x in knots)
    return all(jet(k,k[0],q,"left")==0 and jet(k,k[-1],q,"right")==0 for q in range(3))

def zero_terminal_homogeneous_adjoint_impossible(first_s_weight_rank=3):
    # Z_N=0 and Z_i=Z_{i+1}A_i imply Z_i=0 for every i; hence
    # W_i=Z_{i+1}K_i=0. A regular first S atom is c0(I-K_SS),
    # nonsingular because I-K_SS=R_eff(P_SS+R_eff)^-1.
    return first_s_weight_rank==3

def physical_bounds():
    mean=F(1050297,4840000)
    gamma=F(2432784801508736,10**19)
    floor=F(1,5000)
    return {"sampled_accel_mean_ceiling":mean,
            "physical_joint_vector_floor":gamma,
            "physical_floor_margin":gamma-floor}

def projection_sector_gap():
    # unchanged R_b=.4 and B_a=0.22516660498395405
    return F(2,5)-F(22516660498395405,10**17)

def certificate():
    p=physical_bounds()
    return {
      "qualification":"OU3_SIGNED_TEMPORAL_V1",
      "actual_S_event_signed_temporal_identity":True,
      "homogeneous_zero_terminal_adjoint_closes":False,
      "compatible_adjoint_or_residual_bounds_required":True,
      "sampled_accel_mean_ceiling":str(p["sampled_accel_mean_ceiling"]),
      "physical_joint_vector_floor":str(p["physical_joint_vector_floor"]),
      "physical_floor_margin":str(p["physical_floor_margin"]),
      "projection_sector_gap":str(projection_sector_gap()),
      "regular_A21_pre_projection_BA_precision_ceiling":1000003000,
      "source_uniform_nominal_force_field_temporal_margin":False,
      "source_uniform_nominal_gyro_alias_temporal_margin":False,
      "temporal_margins_imply_finite_B_star":False,
      "six_pivot_floor_implies_finite_B_star":True,
      "carried_adjoint_compatibility_criterion":True,
      "zero_mean_projection_preserves_adjoint":False,
      "signed_physical_bias_summation_by_parts":True,
      "forced_data_adjoint_identity":True,
      "forced_data_adjoint_source_uniform_action_bound":False,
      "physical_span_to_sampled_span":True,
      "construction_uniform_gyro_alias_exclusion":False,
      "B_star_instantiated":False,
      "uniform_historical_AG_readout_action":False,
      "full_21_covariance_upper":False,
      "rho0_certified":False,
    }

if __name__=="__main__":
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))


def separated_reader_action_implication(pivot_floor, coefficient_ceiling,
                                        noise_factor_ceiling, nuisance_ceiling,
                                        terminal_map_ceiling, operation_count):
    """Conditional finiteness from SIX actual residual-norm lower bounds.

    C bounds operator norms of the selected minor, all transports and H_i;
    H bounds the terminal selector and terminal AG map. Both are at least 1.
    p bounds unsquared Gram--Schmidt residual norms; the selector compares
    their squares, which gives the same largest-residual ordering.
    No theorem currently turns Delta_col/Delta_gyr into this pivot premise.
    For six pivots >=p, |det O_I|>=p^6 and ||O_I^-1||<=C^5/p^6.
    The backward recursion starts at the terminal selector, not at zero:
    Y<-Y F or Y<-Y-L_i H_i. Thus ||Y||<=C^N(H+N||L||C).
    Bound complete process/observation factors and the nuisance root only
    AFTER exact AG root cancellation. This coarse scalar ceiling proves
    finiteness; it is not the matrix process comparison used to certify rho.
    """
    from operator import index
    vals=tuple(F(x) for x in (pivot_floor,coefficient_ceiling,noise_factor_ceiling,
                             nuisance_ceiling,terminal_map_ceiling))
    pivot_floor,coefficient_ceiling,noise_factor_ceiling,nuisance_ceiling,terminal_map_ceiling=vals
    operation_count=index(operation_count)
    if (any(x<=0 for x in vals) or coefficient_ceiling<1
            or terminal_map_ceiling<1 or operation_count<1):
        raise ValueError("strict positive source bounds required")
    inv_minor=coefficient_ceiling**5/pivot_floor**6
    reader=terminal_map_ceiling*inv_minor
    y=coefficient_ceiling**operation_count*(terminal_map_ceiling+
                                          operation_count*reader*coefficient_ceiling)
    b_star=(operation_count*noise_factor_ceiling**2*(y**2+reader**2)
            +y**2*nuisance_ceiling)
    return {"pivot_floor":pivot_floor,"inverse_minor_norm_ceiling":inv_minor,
            "reader_norm_ceiling":reader,"B_star_scalar_ceiling":b_star,
            "backward_map_norm_ceiling":y,
            "B_star_finite":True}


def margin_to_Bstar_theorem():
    """Machine-readable status of the controlling implication."""
    return {
      "premises":["inf_W Delta_col(W)>0","inf_W Delta_gyr(W)>0",
                  "shipping coefficient/factor compactness on the fixed finite window"],
      "conclusion":"exists finite B_* with B_W <= B_* I6 on every carried window",
      "rank_structure":"successive largest-residual pivots; observation blocks have rank <=3",
      "proof_gap":"no quantitative six-column pivot bound from the two proposed temporal margins; coefficient compactness is also unproved",
      "implication_closed":False,
      "six_pivot_floor_conditional_implication_closed":True,
      "premise_margins_source_uniformly_certified":False,
    }


def source_margin_attempt():
    """Execute the source-uniform margin arithmetic currently justified.

    The physical part is rigorous.  The nominal transfer remainder is not
    assigned a guessed bound: the shipping assumptions constrain physical
    motion/bias/residuals and recurring magnetic root information, but do not
    yet provide a source-uniform norm ceiling for the signed forced-adjoint
    innovation functional generated by the adaptive gains/resets.  Returning
    None is therefore a mathematical failure location, not infrastructure.
    """
    p=physical_bounds()
    return {
      "physical_collinearity_reserve":p["physical_floor_margin"],
      "physical_gyro_sampling_reserve":F(1,1),  # normalized MAG service premise
      "forced_adjoint_signed_transfer_ceiling":None,
      "reset_transport_transfer_ceiling":None,
      "Delta_col_lower":None,
      "Delta_gyr_lower":None,
      "failure_inequality":"physical reserve - signed forced-adjoint/reset transfer > 0",
      "classification":"missing source-uniform signed transfer bound",
      "new_physical_assumption_needed":False,
    }


def literal_signed_functional_bound(weight_l1, endpoint_state_bound,
                                    affine_defect_l1, affine_defect_bound):
    """Source-uniform bound obtained by exact affine summation by parts.

    For W_i=Z_{i+1}K_i and Z_i=Z_{i+1}A_i the innovation functional is
      sum W_i r_i = Z_N u_N-Z_0 u_0-sum Z_{i+1}d_i.
    CONDITIONAL on both compatibility equations and the stated endpoint and
    defect bounds. Without compatibility the two residual sums remain.
    """
    vals=(weight_l1,endpoint_state_bound,affine_defect_l1,affine_defect_bound)
    if any(x<0 for x in vals): raise ValueError("nonnegative bounds required")
    return 2*weight_l1*endpoint_state_bound + affine_defect_l1*affine_defect_bound


def source_uniform_nominal_endpoint_bounds():
    """Bounds already supplied by literal shipping projection/clamps.

    BA is globally projected.  BG has no analogous projection. AW has a
    covariance/tuner clamp but its *mean* has no shipping saturation. Thus the
    presently proved lemmas do not provide an absolute source-uniform endpoint
    bound for u=(b_hat_g,a_hat_w).  This distinction is the exact reason the
    endpoint telescoping cannot yet become a numeric Delta margin.
    """
    return {
      "b_hat_a_norm":F(2,5),
      "b_hat_g_norm":None,
      "a_hat_w_norm":None,
      "physical_b_g_norm":F(1,50),
      "physical_a_norm":F(44,5),
      "classification":"BG/AW estimator means have no source-uniform absolute clamp",
    }


def forced_adjoint_source_bound():
    """Derive the strongest bound available from the literal same-history recursion."""
    ep=source_uniform_nominal_endpoint_bounds()
    return {
      "identity":"sum W_i r_i = Z_N u_N-Z_0 u_0 + sum (Z_i-Z_(i+1)A_i)u_i + sum (W_i-Z_(i+1)K_i)r_i - sum Z_(i+1)d_i",
      "endpoint_only_identity_requires_compatibility":True,
      "compatibility_source_uniformly_verified":False,
      "innovation_energy_needed":False,
      "independent_gain_box_needed":False,
      "BA_endpoint_bounded":True,
      "BG_endpoint_bounded":ep["b_hat_g_norm"] is not None,
      "AW_endpoint_bounded":ep["a_hat_w_norm"] is not None,
      "finite_numeric_ceiling":None,
      "forced_data_identity_available":True,
      "reason":"homogeneous compatibility residuals can be rewritten by the forced data adjoint, but its root, joint rotation/reference action and literal defects remain unbounded source-uniformly",
      "consequence":"no finite numeric source-uniform ceiling follows from the currently proved lemmas; insufficiency of the physical assumptions is not proved",
    }


def endpoint_annihilating_multiplier_constraints():
    """Candidate boundary conditions, not a carried estimator cancellation.

    The spline jets remove v,p,S boundary terms in integration by parts.
    A continuum gyro companion z_b'=-z_theta has zero end values iff the
    multiplier integral is zero. Neither fact supplies observation forcing
    compatibility or eliminates the discrete residual sums automatically.
    """
    return {
      "AW_endpoint_conditions":["psi(t0)=psi(t1)=0","psi'(t0)=psi'(t1)=0",
                                "psi''(t0)=psi''(t1)=0"],
      "four_S_spline_satisfies_boundary_jets":True,
      "coupled_estimator_endpoint_cancellation_certified":False,
      "BG_companion_equation":"z_b'=-z_theta",
      "BG_endpoint_conditions":["z_b(t0)=0","z_b(t1)=0"],
      "equivalent_BG_moment_condition":"integral z_theta dt = 0",
      "physical_BG_increment_supply":"sum z_b,k w_g,k; ||w_g,k||<=D_g dt_k",
    }


def gyro_zero_mean_companion(z_theta_integrals):
    """Exact discrete companion test for the signed gyro recurrence.

    Inputs are candidate exact cell integrals of an attitude multiplier.
    z_b starts and ends at zero iff their signed sum is zero.  This removes
    a constant component of a physical bias sequence. Compatibility with the
    actual adjoint/observation equations is a separate obligation; it does not
    remove the estimator's constant gyro-bias offset or its alias risk.
    """
    q=tuple(F(x) for x in z_theta_integrals)
    zb=F(0); path=[zb]
    for x in q:
        zb-=x; path.append(zb)
    return {"path":path,"endpoint_zero":zb==0,"zero_mean":sum(q)==0}


def balanced_gyro_weights(cell_integrals):
    """Project one scalar multiplier sequence onto the zero-mean subspace.

    This is an algebraic construction, not a shipping certificate.  It shows
    only that the signed sum vanishes. It generally destroys the carried
    adjoint equations and cannot certify estimator endpoint cancellation.
    """
    q=[F(x) for x in cell_integrals]
    if not q: raise ValueError("cells required")
    mean=sum(q)/len(q)
    b=[x-mean for x in q]
    assert sum(b)==0
    return b


def endpoint_cancelled_source_bound(z_norm_l1, defect_bound,
                                    gyro_companion_l1, gyro_bias_rate,
                                    duration, *, state_residual_bound,
                                    innovation_residual_bound):
    """Conditional bound retaining BOTH compatibility residual supplies.

    Endpoint cancellation, each supplied norm ceiling and its same-history
    interpretation must be proved separately. Neither residual defaults to 0.
    """
    vals=tuple(F(x) for x in (z_norm_l1,defect_bound,gyro_companion_l1,
                              gyro_bias_rate,duration,state_residual_bound,
                              innovation_residual_bound))
    if any(x<0 for x in vals): raise ValueError("nonnegative bounds required")
    z,d,zb,dg,T,rs,ri=vals
    return z*d + zb*dg*T + rs + ri


def adjoint_compatibility(carried_gains, weights):
    """Exact test for Z_N C=W, C_i=A_(N-1)...A_(i+1)K_i.

    Returns a compatible terminal multiplier, or an exact witness v with
    C v=0 and W v!=0. No rank tolerance or normal equations are used.
    All exported-float operands must be converted to their exact rationals;
    such a test then concerns that exported word, not every source history.
    """
    from .matrix_certificates import matmul
    c=[[F(x) for x in row] for row in carried_gains]
    w=[[F(x) for x in row] for row in weights]
    if not c or not w or not c[0] or any(len(r)!=len(c[0]) for r in c+w):
        raise ValueError("nonempty C and W with equal column count required")
    n,m,p=len(c),len(c[0]),len(w)
    basis=[]
    for j in range(m):
        x=[c[i][j] for i in range(n)]
        y=[w[i][j] for i in range(p)]
        v=[F(i==j) for i in range(m)]
        for pivot,b,d,q in basis:
            a=x[pivot]
            if not a: continue
            x=[s-a*t for s,t in zip(x,b)]
            y=[s-a*t for s,t in zip(y,d)]
            v=[s-a*t for s,t in zip(v,q)]
        pivot=next((i for i,a in enumerate(x) if a),None)
        if pivot is None:
            if any(y):
                assert all(a==[0] for a in matmul(c,[[a] for a in v]))
                assert matmul(w,[[a] for a in v])==[[a] for a in y]
                return {"compatible":False,"kernel_witness":v,
                        "weight_residual":y,"carried_rank_at_failure":len(basis)}
        else:
            a=x[pivot]
            basis.append((pivot,[s/a for s in x],[s/a for s in y],[s/a for s in v]))
    z=[[F(0) for _ in range(n)] for _ in range(p)]
    for pivot,b,d,_ in reversed(basis):
        for i in range(p):
            z[i][pivot]=d[i]-sum(z[i][k]*b[k] for k in range(n) if k!=pivot)
    assert matmul(z,c)==w
    return {"compatible":True,"terminal_multiplier":z,"carried_rank":len(basis)}


def signed_bias_terms(weights, biases):
    """Signed Abel identity for one scalar physical bias history.

    Matrix weights follow componentwise. Sum c_i b_i=(sum c_i)b_0
    +sum_j (sum_(i>j)c_i)(b_(j+1)-b_j). Take norms only after these sums.
    """
    c,b=tuple(map(F,weights)),tuple(map(F,biases))
    if not c or len(c)!=len(b): raise ValueError("matching nonempty history required")
    tails=[sum(c[j+1:],F(0)) for j in range(len(c)-1)]
    return {"boundary":sum(c)*b[0],
            "increments":sum((a*(b[j+1]-b[j]) for j,a in enumerate(tails)),F(0)),
            "tail_weights":tails}


def forced_data_adjoint(operations, weights):
    """Exact residual-retaining rewrite, not homogeneous compatibility.

    On ONE actual word, r_i=y_i-H_i u_i+epsilon_i and
    u_(i+1)=A_i u_i+K_i r_i+d_i. For Z_N=0 set
      L_i=W_i+Z_(i+1)K_i,
      Z_i=Z_(i+1)A_i-L_i H_i.
    Then sum W_i r_i=Z_0 u_0+sum L_i(y_i+epsilon_i)
                       +sum Z_(i+1)d_i.
    H is the literal mean observation map, NOT the full EKF Jacobian.
    Rotation/reference dependence, root action and defects remain. This
    enters the tail inequality only through the still-open signed margin.
    """
    from .matrix_certificates import add, matmul
    if not operations or len(operations)!=len(weights):
        raise ValueError("matching nonempty operations and weights required")
    n,p=len(operations[0]['A']),len(weights[0])
    if not n or not p: raise ValueError("nonempty state and output required")
    def cast(a,rows,cols):
        if len(a)!=rows or any(len(row)!=cols for row in a):
            raise ValueError("invalid forced adjoint matrix shape")
        return [[F(x) for x in row] for row in a]
    z=[[F(0) for _ in range(n)] for _ in range(p)]
    states=[z]; inputs=[]
    for op,w in reversed(list(zip(operations,weights))):
        m=len(op['H'])
        if not m: raise ValueError("nonempty innovation coordinates required")
        a=cast(op['A'],n,n); k=cast(op['K'],n,m)
        h=cast(op['H'],m,n); w=cast(w,p,m)
        ell=add(w,matmul(z,k))
        z=add(matmul(z,a),matmul(ell,h),F(-1))
        states.append(z); inputs.append(ell)
    return {"Z":list(reversed(states)),"L":list(reversed(inputs))}


def sampled_tilt_span_lower(theta_e, omega_max, fill_distance):
    """Angular span >= theta_E-2 Omega eta on samples covering the window.

    eta is the physical-time fill distance, not automatically half a step:
    use h_max for in-window samples unless endpoint coverage proves better.
    This is a physical-direction statement, never an estimator-attitude one.
    """
    theta,omega,eta=map(F,(theta_e,omega_max,fill_distance))
    if theta<=0 or omega<0 or eta<0: raise ValueError("invalid span/coverage bounds")
    return max(F(0),theta-2*omega*eta)


def signed_acceleration_supply(endpoint_weight_norm_sum, weight_variation,
                               weighted_cell_square_sum, velocity_bound, jerk_bound):
    """Conditional physical bound on sum h_i C_i a(t_i).

    C_i can be the actual signed accelerometer residual multiplier times the
    true body-to-world inverse. First sum by parts against physical velocity:
    C_last v_N-C_first v_0+sum (C_(i-1)-C_i)v_i. Left-cell sampling adds
    at most (J/2)sum ||C_i||h_i^2. No estimator mean or innovation is bounded
    by this lemma. Norms are taken after signed temporal weight differences.
    """
    e,var,cells,v,j=map(F,(endpoint_weight_norm_sum,weight_variation,
                          weighted_cell_square_sum,velocity_bound,jerk_bound))
    if min(e,var,cells,v,j)<0: raise ValueError("nonnegative source bounds required")
    return v*(e+var)+j*cells/2


def physical_interval_sampling_error(weight_norms, times, left, right, jerk_bound):
    """Jerk remainder after signed weights are summed within one interval.

    D=sum C_i, H=right-left. The exact physical relation is
      sum C_i a(t_i)=(D/H)(v(right)-v(left))+epsilon.
    Integrating ||a(t_i)-a(s)||<=J|t_i-s| bounds epsilon by
      J sum ||C_i|| ((t_i-left)^2+(right-t_i)^2)/(2H).
    D MUST be summed with signs before its operator norm is bounded. The
    intervals use one carried physical velocity history; no restart or
    synthetic estimator update is involved.
    """
    weights=tuple(map(F,weight_norms)); times=tuple(map(F,times))
    left,right,jerk=map(F,(left,right,jerk_bound))
    if (not weights or len(weights)!=len(times) or left>=right or jerk<0
            or any(w<0 for w in weights) or any(t<left or t>right for t in times)):
        raise ValueError("valid physical interval, samples and norms required")
    return jerk*sum((w*((t-left)**2+(right-t)**2)/(2*(right-left))
                     for w,t in zip(weights,times)),F(0))


def gyro_construction_barrier():
    """Necessary bias-mean magnitude for a first complete turn; not exclusion."""
    import json
    from pathlib import Path
    c=json.loads(Path(__file__).with_name('constants.json').read_text(),parse_float=F)
    omega=c['marine_motion']['Omega_max_rad_s']; bg=c['imu_bias']['B_g_rad_s']
    ng=c['sensor_model']['gyro_fast_residual_norm_max_rad_s']
    hmax=c['sensor_model']['sample_period_max_s']; pi_lower=F('3.14159265358979323846')
    return {"first_prediction_increment_ceiling":hmax*(omega+bg+ng),
            "complete_turn_requires_bias_norm_at_least":2*pi_lower/hmax-omega-bg-ng,
            "initial_bias_mean":F(0),
            "all_time_signed_gain_innovation_sum_ceiling":None,
            "construction_unreachable_certified":False}
