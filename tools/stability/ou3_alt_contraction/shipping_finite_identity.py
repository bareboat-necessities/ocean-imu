#!/usr/bin/env python3
"""Build passive shipping probes and verify finite-word correspondence.

This is an implementation REGRESSION, not a source-uniform certificate or a
storage-feasibility experiment.  The executable runs the actual wrapper and
its private frontend.  A temporary include overlay adds observation calls only;
the tracked shipping headers are never rewritten.  A second, uninstrumented
build runs the same input and its recorded sample states must match byte-for-
byte.  The probe does not set filter states, gains, schedules, or mode guards.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[3]
CORE = Path('src/kalman_ou_iii/Kalman3D_Wave_OU_III.h')
WRAPPER = Path('src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h')


def function_span(text: str, signature: str) -> tuple[int, int]:
    """Find a unique function body, respecting comments and quoted strings."""
    if text.count(signature) != 1:
        raise ValueError(f'expected one function signature: {signature!r}')
    start = text.index('{', text.index(signature))
    depth, i, state = 1, start+1, 'code'
    while i < len(text):
        ch, nxt = text[i], text[i:i+2]
        if state == 'line':
            if ch == '\n': state = 'code'
        elif state == 'block':
            if nxt == '*/': state, i = 'code', i+1
        elif state in ('"', "'"):
            if ch == '\\': i += 1
            elif ch == state: state = 'code'
        elif nxt == '//': state, i = 'line', i+1
        elif nxt == '/*': state, i = 'block', i+1
        elif ch in ('"', "'"): state = ch
        elif ch == '{': depth += 1
        elif ch == '}':
            depth -= 1
            if not depth: return start, i
        i += 1
    raise ValueError(f'unclosed function body: {signature}')


def erase_observation_calls(text: str) -> str:
    # Delete only the inserted call, not an original trailing source comment.
    pattern = (r'ou3_alt_probe::(?:record|innovation)\([^;\n]*\);|'
               r'ou3_alt_probe::Scope<decltype\(\*this\)> alt_scope\([^;\n]*\);')
    return re.sub(pattern, '', text)


def instrument_core(text: str) -> tuple[str, list[dict]]:
    """Insert passive calls at literal operation boundaries, fail on drift."""
    if 'ou3_alt_probe::' in text:
        raise ValueError('source already contains probe calls')
    original = text
    manifest = []

    def edit(sig, needle=None, insert='', *, at_entry=False):
        nonlocal text
        lo, hi = function_span(text, sig)
        body = text[lo:hi+1]
        if at_entry:
            updated = body[:1]+'\n'+insert+'\n'+body[1:]
        else:
            if body.count(needle) != 1:
                raise ValueError(f'nonunique/missing hook {sig}: {needle!r}')
            updated = body.replace(needle, needle+'\n'+insert, 1)
        manifest.append({'signature': sig, 'anchor': 'entry' if at_entry else needle,
                         'insertion': insert})
        text = text[:lo]+updated+text[hi+1:]

    tu = '::time_update('\
        '\n    Vector3 const& gyr_body, T Ts)'
    edit(tu, at_entry=True, insert='ou3_alt_probe::record(1, *this);')
    edit(tu, '    apply_pending_aw_covariance_inflation_();',
         '    ou3_alt_probe::record(3, *this);')
    # Prediction must be recorded BEFORE the floor, not after it.
    lo, hi = function_span(text, tu)
    body = text[lo:hi+1]
    anchor = '    apply_pending_aw_covariance_inflation_();'
    body = body.replace(anchor, '    ou3_alt_probe::record(2, *this);\n'+anchor)
    text = text[:lo]+body+text[hi+1:]
    manifest.append({'signature': tu, 'anchor': 'before '+anchor,
                     'insertion': 'ou3_alt_probe::record(2, *this);'})
    edit(tu, '    symmetrize_Pext_();', '    ou3_alt_probe::record(4, *this);')

    for kind, sig in ((10, '::measurement_update_acc_only('),
                      (20, '::measurement_update_mag_only('),
                      (30, '::applyIntegralZeroPseudoMeas()')):
        edit(sig, at_entry=True,
             insert=f'ou3_alt_probe::Scope<decltype(*this)> alt_scope({kind}, *this);')
        lo, hi = function_span(text, sig)
        body = text[lo:hi+1]
        anchor = 'xext.noalias() += K * r;'
        if body.count(anchor) != 1:
            raise ValueError('mean-update hook lost')
        body = body.replace(anchor,
            f'ou3_alt_probe::innovation({kind+2}, *this, r, PCt, S_mat, K, '+
            ('acc_meas' if kind == 10 else 'mag_meas' if kind == 20 else 'Vector3::Zero()')+');\n    '+
            anchor+f'\n    ou3_alt_probe::record({kind+3}, *this);', 1)
        anchor2 = 'joseph_update3_(K, S_mat, PCt);'
        if body.count(anchor2) != 1:
            raise ValueError('Joseph hook lost')
        body = body.replace(anchor2, anchor2+f'\n    ou3_alt_probe::record({kind+4}, *this);', 1)
        text = text[:lo]+body+text[hi+1:]
        manifest.append({'signature': sig, 'anchor': 'mean and Joseph',
                         'insertion': 'innovation / post-mean / post-Joseph observation'})

    reset = '::applyQuaternionCorrectionFromErrorState()'
    edit(reset, at_entry=True, insert='ou3_alt_probe::Scope<decltype(*this)> alt_scope(40, *this);')
    edit(reset, '    qref.normalize();', '    ou3_alt_probe::record(42, *this);')
    edit(reset, '    apply_error_state_reset_jacobian_(dtheta);',
         '    ou3_alt_probe::record(43, *this);')
    edit(reset, '    // Clear the local attitude-error state after injection.\n    xext.template head<3>().setZero();',
         '    ou3_alt_probe::record(44, *this);')
    edit('::project_acc_bias_()', at_entry=True,
         insert='ou3_alt_probe::Scope<decltype(*this)> alt_scope(50, *this);')
    edit('void set_acc_bias_updates_enabled(bool en)', at_entry=True,
         insert='ou3_alt_probe::Scope<decltype(*this)> alt_scope(60, *this);')
    erased = erase_observation_calls(text)
    if ''.join(erased.split()) != ''.join(original.split()):
        raise AssertionError('observation erasure did not recover shipping source')
    return text, manifest


def make_overlay(directory: Path) -> dict:
    original = (ROOT/CORE).read_text()
    transformed, patches = instrument_core(original)
    target = directory/'kalman_ou_iii'
    target.mkdir(parents=True, exist_ok=True)
    (target/CORE.name).write_text(transformed)
    shutil.copyfile(ROOT/WRAPPER, target/WRAPPER.name)
    return {'core_sha256': hashlib.sha256(original.encode()).hexdigest(),
            'wrapper_sha256': hashlib.sha256((ROOT/WRAPPER).read_bytes()).hexdigest(),
            'overlay_sha256': hashlib.sha256(transformed.encode()).hexdigest(),
            'patches': patches}


def eigen_include() -> Path:
    candidates = [os.environ.get('EIGEN_INCLUDE_DIR', ''), '/usr/include/eigen3',
                  '/usr/local/include/eigen3',
                  '/opt/pyvenv/lib/python3.13/site-packages/casadi/include/eigen3']
    for s in candidates:
        p = Path(s)
        if s and (p/'Eigen/Core').is_file(): return p
    raise RuntimeError('Eigen3 not found; set EIGEN_INCLUDE_DIR')


def same_files(a: Path, b: Path) -> bool:
    if a.stat().st_size != b.stat().st_size: return False
    with a.open('rb') as fa, b.open('rb') as fb:
        while True:
            aa, bb = fa.read(1 << 20), fb.read(1 << 20)
            if aa != bb: return False
            if not aa: return True


def run(output: Path, work: Path) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    overlay = work/'include'
    manifest = make_overlay(overlay)
    source_before = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                     for p in (CORE, WRAPPER)}
    compiler = os.environ.get('CXX', 'g++')
    common = [compiler, '-std=c++20', '-O1', '-fno-fast-math', '-ffp-contract=off',
              '-DEIGEN_NON_ARDUINO', '-I'+str(eigen_include())]
    source = ROOT/'tests/ou3_alt_contraction/shipping_finite_word.cpp'
    # Link only the canonical gravity definition; the unrelated simulation
    # translation unit is not part of this host correspondence harness.
    gravity_source = ROOT/'src/util/W3dSimCommon.cpp'
    definitions = re.findall(r'^extern const float g_std = [^;]+;',
                             gravity_source.read_text(), re.MULTILINE)
    if len(definitions) != 1:
        raise RuntimeError('canonical g_std definition changed')
    gravity_unit = work/'gravity.cpp'
    gravity_unit.write_text(definitions[0]+'\n')
    manifest['gravity_definition'] = definitions[0]
    manifest['gravity_source_sha256'] = hashlib.sha256(gravity_source.read_bytes()).hexdigest()
    timings, reports = {}, {}
    for label, observed in (('plain', False), ('observed', True)):
        exe = work/label
        command = common + (['-I'+str(overlay)] if observed else []) + [
            '-I'+str(ROOT/'src'), '-DALT_OBSERVE='+str(int(observed)),
            str(source), str(gravity_unit), '-o', str(exe)]
        start = time.perf_counter()
        subprocess.run(command, check=True, cwd=ROOT)
        timings[label+'_compile_s'] = time.perf_counter()-start
        dest = work/(label+'-data')
        dest.mkdir(exist_ok=True)
        start = time.perf_counter()
        rss = work/(label+'-rss-kib.txt')
        command = [str(exe), str(dest)]
        if Path('/usr/bin/time').is_file():
            command = ['/usr/bin/time', '-f', '%M', '-o', str(rss), *command]
        proc = subprocess.run(command, check=True, capture_output=True, text=True)
        timings[label+'_peak_rss_kib'] = int(rss.read_text().strip()) if rss.is_file() else None
        timings[label+'_run_s'] = time.perf_counter()-start
        reports[label] = json.loads(proc.stdout)
    unchanged = same_files(work/'plain-data/samples.bin', work/'observed-data/samples.bin')
    if not unchanged:
        raise AssertionError('passive instrumentation changed recorded shipping sample state')
    if any(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != d
           for p, d in source_before.items()):
        raise AssertionError('tracked shipping source changed')
    # Imported late: the build/trace audit remains usable independently.
    from tools.stability.ou3_alt_contraction.verify_finite_trace import verify
    verification = verify(work/'observed-data')
    result = {'qualification': 'FINITE_SHIPPING_WORD_IMPLEMENTATION_REGRESSION_ONLY',
              'source': manifest, 'timings': timings, 'runs': reports,
              'recorded_sample_states_bit_identical': unchanged,
              'tracked_shipping_headers_unchanged': True,
              'observation_calls_erase_to_original_source': True,
              'frontend_audit_is_not_complete_state_or_source_qualification': True,
              'finite_identity_regression': verification,
              'source_uniform_word_qualified': False,
              'storage_feasibility_attempted': False,
              'deployment_roundoff_enclosed': False,
              'ALT_LIVE_PASS': False, 'ALT_STARTUP_PASS': False,
              'ALT_END_TO_END_PASS': False}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path)
    args = parser.parse_args()
    if args.work_dir:
        result = run(args.output, args.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='ou3-alt-finite-') as tmp:
            result = run(args.output, Path(tmp))
    print(json.dumps({k: result[k] for k in ('qualification',
                     'recorded_sample_states_bit_identical',
                     'finite_identity_regression', 'timings')}, indent=2))
    return 0


if __name__ == '__main__':
    # Running by pathname must not depend on the caller's PYTHONPATH.
    import sys
    sys.path[:0] = [str(ROOT), str(ROOT/'tools/stability')]
    raise SystemExit(main())
