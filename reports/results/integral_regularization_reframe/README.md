# Integral-regularization paper reframe

This branch intentionally changes publication claims, not production filter tuning.

The standalone paper now treats OU similarity and the scalar drift-band model as scaling/mechanism arguments. The reduced drift--distortion minimizer is explicitly illustrative and is not presented as the optimum of the complete attitude/bias-coupled MEKF. The deployed schedule remains empirically justified.

Horizontal conclusions are reduced to the supported claim that isotropic integral regularization performs adequately on the tested directional records. Separate world-axis optima are not claimed because surge/sway projection depends on wave direction and the available heading sample is limited.

The publication contract in `tests/validation/test_integral_regularization_paper.py` prevents reintroduction of the removed optimum, legacy horizontal-tuning, and table claims.
