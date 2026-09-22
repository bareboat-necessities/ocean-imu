# Independent horizontal R_S split grid (training seeds)

30 arms x 8 pinned vessel-RAO records x 3 training seeds = 720 paired cases,
c_tau and c_sigma held at the pre-RAO operating point, Q_bg at 3.75e-10.
Ratios are the arm over the paired baseline case; below 1 is better.

## The two horizontal axes control different attitude channels

Pitch RMS ratio, by k_Sy (columns) and k_Sx (rows):

|k_Sx\k_Sy| 0.40 | 0.50 | 0.65 | 0.85 | 1.15 | 1.45 |
|---|---|---|---|---|---|---|
|0.50|1.714|1.272|1.093|1.062|1.038|1.027|
|0.72|1.810|1.355|1.098|1.021|1.001|1.011|
|0.90|1.841|1.386|1.111|1.021|0.999|1.017|
|1.15|1.861|1.407|1.122|1.023|1.000|1.019|
|1.45|1.873|1.417|1.128|1.022|0.997|1.018|

Roll RMS ratio, same layout:

|k_Sx\k_Sy| 0.40 | 0.50 | 0.65 | 0.85 | 1.15 | 1.45 |
|---|---|---|---|---|---|---|
|0.50|1.041|1.014|0.994|0.982|0.977|0.979|
|0.72|1.051|1.026|1.006|0.995|0.991|0.992|
|0.90|1.056|1.031|1.011|1.000|0.996|0.997|
|1.15|1.060|1.035|1.015|1.004|1.000|1.001|
|1.45|1.063|1.038|1.019|1.008|1.003|1.004|

Pitch is set by k_Sy and is nearly flat in k_Sx; roll is set by k_Sx and is
nearly flat in k_Sy. Pitch is the steep one: it degrades 71 to 87 percent at
k_Sy = 0.40 and reaches its minimum at k_Sy = 1.15, the value TFG already
ships. Roll moves less than 6 percent across the whole grid.

## Deployed OU-III's split is strongly wrong for TFG

OU-III ships k_Sx = 0.72, k_Sy = 0.50 on these same records, tightening the
keel-damped sway axis. Applied to TFG that arm ranks 22 of 30: pitch +35.5
percent, accelerometer bias +9.5 percent, 3-D displacement +7.0 percent,
physical score +0.059. TFG wants the opposite anisotropy, k_Sx < k_Sy, and the
low-k_Sy corner is the worst region of the entire grid.

The keel reasoning that justifies OU-III's split is sound about the dataset,
but it does not transfer: the two filters are not failing on the same axis.

## What the split is actually worth

The best arms are interior in both axes, not on a face:

| arm | k_Sx | k_Sy | physical score | pitch | accel bias | 3-D | gate-fail cases |
|---|---|---|---|---|---|---|---|
| split_x0.72_y1.15 | 0.72 | 1.15 | -0.00386 | +0.07% | -1.15% | -0.73% | 19 |
| split_x0.9_y1.15 | 0.90 | 1.15 | -0.00328 | -0.05% | -0.48% | -0.80% | 19 |
| split_x0.72_y1.45 | 0.72 | 1.45 | -0.00100 | +1.10% | -1.30% | -0.20% | 20 |
| baseline | 1.15 | 1.15 | 0.00000 | - | - | - | 16 |

k_Sx = 0.50 is worse than 0.72, so the x optimum is interior near 0.72 to 0.90.
k_Sy = 1.45 is worse than 1.15, so the y optimum is interior at about 1.15.

This independently reproduces the extended four-dimensional grid's report of an
interior x optimum near 0.75 with y above it, from a separate execution over a
much wider region, and it locates the y optimum at the shipping value rather
than at 1.3 to 1.4.

The gain does not survive its own consistency check. The leading arm improves
only 15 of 24 paired cases, so the 0.4 percent composite is not a consistent
effect across records and seeds, and every anisotropic arm adds gate-failure
cases against the baseline's 16. Nothing here promotes.
