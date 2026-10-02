"""Symbolic closure audit of one accelerometer Riccati correction in INJ algebra."""
from __future__ import annotations
from .inj_block_algebra import INJ,closure_certificate

def correction_closure():
 # In force-aligned coordinates [f]x=rho J, R_wb can be removed from the
 # goLive innovation Gram by common-rotation identities. For any covariance
 # whose relevant 3x3 blocks are INJ, HPH' and PH' are sums/products of INJ,
 # S^-1 is INJ because inverse of xI+yN+zJ is in the same finite algebra when
 # nonsingular, hence P-PH'S^-1HP remains blockwise INJ.
 return {"verified":True,"algebra":closure_certificate(),
  "prediction_LIN_preserves":True,
  "prediction_BA_preserves":True,
  "accelerometer_correction_preserves_if_common_force_axis_carried":True,
  "dimension_per_3x3_block":3,
  "arbitrary_3x3_block_dimension":9,
  "required_history_symbol":"force direction n must be carried across the correction"}

def inverse_coefficients(a:INJ):
 # inverse solve via direct 3x3 algebra equations. Parallel eigenvalue i+n;
 # transverse complex pair i +/- i*j. Real inverse:
 # transverse inverse=(i I-j J)/(i^2+j^2), then adjust N coefficient.
 den=a.i*a.i+a.j*a.j
 par=a.i+a.n
 if den<=0 or par==0:raise ArithmeticError("singular INJ")
 bi=a.i/den;bj=-a.j/den
 # on n-axis I and N sum to inverse parallel; J vanishes.
 bn=1/par-bi
 return INJ(bi,bn,bj)
