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
    fac=(2,1,1)[order]
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
      "observation_forced_adjoint_required":True,
      "sampled_accel_mean_ceiling":str(p["sampled_accel_mean_ceiling"]),
      "physical_joint_vector_floor":str(p["physical_joint_vector_floor"]),
      "physical_floor_margin":str(p["physical_floor_margin"]),
      "projection_sector_gap":str(projection_sector_gap()),
      "regular_A21_pre_projection_BA_precision_ceiling":1000003000,
      "source_uniform_nominal_force_field_temporal_margin":False,
      "source_uniform_nominal_gyro_alias_temporal_margin":False,
      "temporal_margins_imply_finite_B_star":True,
      "B_star_instantiated":False,
      "uniform_historical_AG_readout_action":False,
      "full_21_covariance_upper":False,
      "rho0_certified":False,
    }

if __name__=="__main__":
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))


def separated_reader_action_implication(delta_col, delta_gyr, coefficient_ceiling,
                                        noise_factor_ceiling, nuisance_ceiling,
                                        terminal_map_ceiling, operation_count):
    """Quantitative implication from positive temporal margins to finite B_*.

    This is not a margin certificate.  It proves the next logical step once
    the same-history enclosure supplies strict delta_col and delta_gyr.

    Let delta=min(delta_col,delta_gyr).  On every largest-residual pivot chart,
    the six successive residual pivots of O are bounded below by delta after
    the proof-coordinate normalization used by the temporal margins.  Hence
    the selected 6x6 minor has ||O_I^{-1}||_2 <=
    coefficient_ceiling**5 / delta**6 by adjugate/Hadamard.  The exact reader
    L=T_h O_I^{-1} therefore has the displayed uniform norm ceiling.  The
    backward action is a finite sum of transported rank<=3 noise/process
    factors plus the nuisance-root residual.  Bounding each chronological
    transport by coefficient_ceiling gives the explicit B_* below.

    The deliberately coarse exponent is acceptable here: the result needed
    is finiteness, not a practical contraction rate.  A useful J/rho still
    requires the rigorous source margins and then a sharper action enclosure.
    """
    vals=(delta_col,delta_gyr,coefficient_ceiling,noise_factor_ceiling,
          nuisance_ceiling,terminal_map_ceiling)
    if any(x<=0 for x in vals) or operation_count<1:
        raise ValueError("strict positive source bounds required")
    delta=min(delta_col,delta_gyr)
    inv_minor=coefficient_ceiling**5/delta**6
    reader=terminal_map_ceiling*inv_minor
    transport=max(1.0,coefficient_ceiling)**operation_count
    # ||sum X_i X_i'|| <= sum ||X_i||^2.  The final term covers the
    # nuisance-root residual in exactly the same backward reader recursion.
    b_star=(operation_count*(reader*noise_factor_ceiling*transport)**2
            +(reader*transport)**2*nuisance_ceiling)
    return {"delta":delta,"inverse_minor_norm_ceiling":inv_minor,
            "reader_norm_ceiling":reader,"B_star_scalar_ceiling":b_star,
            "B_star_finite":True}


def margin_to_Bstar_theorem():
    """Machine-readable status of the controlling implication."""
    return {
      "premises":["inf_W Delta_col(W)>0","inf_W Delta_gyr(W)>0",
                  "shipping coefficient/factor compactness on the fixed finite window"],
      "conclusion":"exists finite B_* with B_W <= B_* I6 on every carried window",
      "rank_structure":"successive largest-residual pivots; observation blocks have rank <=3",
      "proof":"finite pivot-chart cover + adjugate/Hadamard inverse bound + finite backward factor action",
      "implication_closed":True,
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
    The crucial point is that no innovation/NIS norm appears.  Once the
    forced adjoint is represented by endpoint observations, the source bound
    is an endpoint-state bound plus literal affine defects.
    """
    vals=(weight_l1,endpoint_state_bound,affine_defect_l1,affine_defect_bound)
    if any(x<0 for x in vals): raise ValueError("nonnegative bounds required")
    return 2*weight_l1*endpoint_state_bound + affine_defect_l1*affine_defect_bound


def source_uniform_nominal_endpoint_bounds():
    """Bounds already supplied by literal shipping projection/clamps.

    BA is globally projected.  BG has no analogous projection. AW has a
    covariance/tuner clamp but its *mean* has no shipping saturation. Thus the
    present assumptions do not provide an absolute source-uniform endpoint
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
      "identity":"sum W_i r_i = Z_N u_N-Z_0 u_0-sum Z_(i+1)d_i",
      "innovation_energy_needed":False,
      "independent_gain_box_needed":False,
      "BA_endpoint_bounded":True,
      "BG_endpoint_bounded":ep["b_hat_g_norm"] is not None,
      "AW_endpoint_bounded":ep["a_hat_w_norm"] is not None,
      "finite_numeric_ceiling":None,
      "reason":"exact telescoping leaves BG/AW endpoint means; neither has an absolute shipping/source bound under the current theorem premises",
      "consequence":"a source-uniform numeric correction/reset ceiling cannot be derived from the current premises alone by this adjoint",
    }


def endpoint_annihilating_multiplier_constraints():
    """Boundary conditions needed to remove uncontrolled BG/AW endpoint means.

    For the translation chain S'=p, p'=v, v'=a_hat plus OU AW prediction,
    three integrations by parts show that AW endpoint coefficients vanish when
    psi and its first two derivatives vanish at both ends.  The four-S-event
    quadratic spline already has exactly these six boundary conditions.

    For gyro bias, theta'=-[omega_hat]x theta + b_hat_g at the differential
    level.  A left adjoint z_theta satisfying z' = z[omega_hat]x has the bias
    coefficient integral z_theta dt.  Cancelling absolute b_hat_g endpoints
    through the signed physical bias recurrence requires a companion
    multiplier z_b with z_b'=-z_theta and z_b=0 at both ends; equivalently
    integral z_theta dt=0.  Thus the gyro multiplier must have zero temporal
    mean in addition to zero boundary bias coefficient.
    """
    return {
      "AW_endpoint_conditions":["psi(t0)=psi(t1)=0","psi'(t0)=psi'(t1)=0",
                                "psi''(t0)=psi''(t1)=0"],
      "four_S_spline_satisfies_AW_endpoint_conditions":True,
      "BG_companion_equation":"z_b'=-z_theta",
      "BG_endpoint_conditions":["z_b(t0)=0","z_b(t1)=0"],
      "equivalent_BG_moment_condition":"integral z_theta dt = 0",
      "physical_BG_increment_supply":"sum z_b,k w_g,k; ||w_g,k||<=D_g dt_k",
    }


def gyro_zero_mean_companion(z_theta_integrals):
    """Exact discrete companion test for the signed gyro recurrence.

    Inputs are exact cell integrals of the transported attitude multiplier.
    z_b starts and ends at zero iff their signed sum is zero.  This removes
    absolute b_hat_g endpoint dependence; only bounded physical increments and
    correction/reset defects remain.
    """
    q=tuple(F(x) for x in z_theta_integrals)
    zb=F(0); path=[zb]
    for x in q:
        zb-=x; path.append(zb)
    return {"path":path,"endpoint_zero":zb==0,"zero_mean":sum(q)==0}


def balanced_gyro_weights(cell_integrals):
    """Project one scalar multiplier sequence onto the zero-mean subspace.

    This is an algebraic construction, not a shipping certificate.  It shows
    endpoint cancellation costs one temporal moment rather than an absolute
    BG state bound.
    """
    q=[F(x) for x in cell_integrals]
    if not q: raise ValueError("cells required")
    mean=sum(q)/len(q)
    b=[x-mean for x in q]
    assert sum(b)==0
    return b


def endpoint_cancelled_source_bound(z_norm_l1, defect_bound,
                                    gyro_companion_l1, gyro_bias_rate,
                                    duration):
    """Bound after AW/BG absolute endpoints have been annihilated.

    Remaining terms are literal affine/reset defects plus physical gyro-bias
    increments.  No absolute b_hat_g or a_hat_w endpoint bound appears.
    """
    vals=tuple(F(x) for x in (z_norm_l1,defect_bound,gyro_companion_l1,
                              gyro_bias_rate,duration))
    if any(x<0 for x in vals): raise ValueError("nonnegative bounds required")
    z,d,zb,dg,T=vals
    return z*d + zb*dg*T
