# Regression evidence

`limits.json` and `per_record_limits.json` are the fixed v1.4.0 limits for both
following and head seas. Current measurements and charts live in
[`../sim_v1.4.0/`](../sim_v1.4.0/). Normal CI reads these limits without fitting them.
The steering-separated native baseline updates every horizontal profile and
adds an independent Y cap for each of its 272 cases. See
[`heading_limit_update.md`](heading_limit_update.md) for the new ceilings,
per-record Y examples and explicitly disclosed bias/attitude cap increases.

The other JSON/CSV bundles in this folder retain their original v1.3.1 inputs,
executables and experimental scope: stationary/underway comparisons, stillness
controls, CPU reproducibility and prior study artifact receipts.
`stationary_comparison_inputs.json` explicitly pins the old moving arm used by
the historical encounter-regression inventory. These are not v1.4.0 scores.
