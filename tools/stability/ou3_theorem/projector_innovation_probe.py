"""Verify projector inverse against literal goLive Racc interval."""
import json
from pathlib import Path
from .projector_innovation_inverse import direction_free_bounds,inverse_coefficients
from .release_interval_propagation import one_sample_boxes
ROOT=Path(__file__).resolve().parents[3]
def diagnostic():
 _,_,_,R,*_=one_sample_boxes()
 # Literal Racc is isotropic diagonal interval.
 rlo=R.mid[0][0]-R.rad[0][0];rhi=R.mid[0][0]+R.rad[0][0]
 z=direction_free_bounds(0.,18.60665,rlo,rhi)
 return {**z,"Racc_variance_lower":rlo,"Racc_variance_upper":rhi,
  "rho_max":18.60665,"rho_max_coefficients":inverse_coefficients(18.60665,rlo)}
