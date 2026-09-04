from pathlib import Path

parity = Path("tests/validation/test_ou_paper_code_parity.py")
text = parity.read_text(encoding="utf-8")
old = '            "const float settle_sec = 6.0f / lambda_;",'
new = '''            "const float moment_start_sec = 3.0f / lambda_;",
            "const float usable_floor_sec = 4.0f / lambda_;",
            "const float settled_floor_sec = 6.0f / lambda_;",
            "bool hasUsablePeriod() const",'''
if text.count(old) != 1:
    raise SystemExit("paper/code parity settle token changed unexpectedly")
parity.write_text(text.replace(old, new), encoding="utf-8")

ledger = Path("docs/ou3-proof-research-state.md")
text = ledger.read_text(encoding="utf-8")
old = '''`tools/ou3_sea3_wave_period_frontend.py` remains replay free and non-promoting.
It source-certifies the fixed 0.2 Hz prior, prior-to-first-finite-estimator
handoff, one-sample tuner/period-estimator ordering edge, and the fact that
filter `TunerReady` is not WavePeriodEstimator readiness. Its single-frequency
discrete front-end certificate keeps period warping below about 59 ppm on the
committed 5 ms / 0.03--1.2 Hz channel.'''
new = '''`tools/ou3_sea3_wave_period_frontend.py` remains replay free and non-promoting.
It source-certifies the fixed 0.2 Hz prior, the same estimator's `3/lambda`
moment-start, `4/lambda` plus one-period startup-usable gate, strict
`6/lambda` readiness floor, and the one-sample tuner/period-estimator ordering
edge. OU-II, OU-III, and TFG now require that startup-usable measured period
before Live entry while keeping `WavePeriodEstimator::isReady()` as the stricter
diagnostic state. Its single-frequency discrete front-end certificate keeps
period warping below about 59 ppm on the committed 5 ms / 0.03--1.2 Hz channel.'''
if text.count(old) != 1:
    raise SystemExit("proof-ledger startup paragraph changed unexpectedly")
ledger.write_text(text.replace(old, new), encoding="utf-8")
