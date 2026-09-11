"""Observation-only attachment to actual shipping operations; no proof promotion.

The temporary header overlay adds callbacks but changes no shipping expression.
Anchor counts fail closed on changed attachment sites; source hashes record provenance. The C++
probe forks the whole process at word entry, retaining all hidden state.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
HEADER = 'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'


def overlay(text: str) -> str:
    """Exact function-local anchors; no changes to the deployed header."""
    def insert(old, new, count=1):
        nonlocal text
        if text.count(old) != count:
            raise ValueError(f'shipping attachment anchor changed: {old[:90]!r}')
        text = text.replace(old, new)
    # Scope destructors also observe early-return/rejected paths.
    for signature, tag in (
        ('Vector3 const& gyr_body, T Ts)', 'time_end'),
        ('Vector3 const& acc_meas_body, T tempC)', 'acc_end'),
        ('const Vector3& mag_meas_body)', 'mag_end'),
        ('void Kalman3D_Wave_OU_III<T, with_gyro_bias, with_accel_bias>::applyIntegralZeroPseudoMeas()', 'S_end'),
    ):
        old=signature+'\n{'
        insert(old, old+f'\n    AltScope alt_scope("{tag}", *this);')
    insert('    apply_pending_aw_covariance_inflation_();\n    symmetrize_Pext_();   // Symmetry hygiene',
           '    alt_state("prediction", *this);\n    apply_pending_aw_covariance_inflation_();\n'
           '    symmetrize_Pext_();   // Symmetry hygiene\n    alt_state("aw_floor", *this);')
    # These three solve graphs include actual masked PCt, actual S and actual K.
    for tag, begin_marker, end_marker in (
        ('acc', '::measurement_update_acc_only(', '::measurement_update_mag_only('),
        ('mag', '::measurement_update_mag_only(', '// specific force prediction'),
        ('S', '::applyIntegralZeroPseudoMeas()', '::measurement_update_position_pseudo('),
    ):
        a=text.index(begin_marker); b=text.index(end_marker,a+len(begin_marker))
        segment=text[a:b]
        anchor='    xext.noalias() += K * r;'
        if segment.count(anchor)!=1: raise ValueError(f'{tag} correction anchor changed')
        line=next(line for line in segment.splitlines(keepends=True) if anchor in line)
        segment=segment.replace(line,f'    alt_solve("{tag}", *this, r);\n'+line+
                                '    alt_solve_finish(*this);\n')
        anchor2='    applyQuaternionCorrectionFromErrorState();'
        if segment.count(anchor2)!=1: raise ValueError(f'{tag} injection anchor changed')
        segment=segment.replace(anchor2,'    alt_state("joseph", *this);\n'+anchor2)
        text=text[:a]+segment+text[b:]
    insert('    project_acc_bias_();\n}',
           '    alt_state("reset", *this);\n    project_acc_bias_();\n    alt_state("projection", *this);\n}')
    return text



def strip_observations(text: str) -> str:
    """An audit check: removing only inserted callback lines recovers shipping."""
    return ''.join(line for line in text.splitlines(keepends=True)
                   if not any(token in line for token in
                              ('    AltScope alt_scope(', '    alt_state(', '    alt_solve(', '    alt_solve_finish(')))


def build(directory: Path, eigen: Path, instrumented: bool=True) -> Path:
    src=ROOT/HEADER
    original=src.read_text()
    target=directory/'include/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
    target.parent.mkdir(parents=True,exist_ok=True)
    overlay_text=overlay(original) if instrumented else original
    if strip_observations(overlay_text) != original:
        raise ValueError("observation overlay changed a shipping expression")
    target.write_text(overlay_text)
    # SeaState's include is src-relative, so shadow it without changing a byte.
    front=target.with_name('SeaStateFusionFilter_OU_III.h')
    front.write_bytes((ROOT/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h').read_bytes())
    binary=directory/'shipping-word-probe'
    command=[os.environ.get('CXX','g++'),'-std=c++20','-O2','-fno-fast-math',
             '-I'+str(directory/'include'),'-I'+str(ROOT/'src'),'-isystem',str(eigen),
             str(ROOT/'tests/ou3_alt_contraction/shipping_word_probe.cpp'),'-o',str(binary)]
    subprocess.run(command,check=True)
    (directory/'provenance.json').write_text(json.dumps({
        'shipping_header_sha256':hashlib.sha256(original.encode()).hexdigest(),
        'overlay_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
        'compile_command':command, 'shipping_sources_modified':False,'instrumented':instrumented,
        'stripping_callbacks_recovers_original_bytes':True,
        'scope':'actual binary32 paired execution; not source-uniform qualification',
        'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    },indent=2)+'\n')
    return binary


def run(directory: Path, binary: Path, family: int, epsilon: float, detailed: int):
    directory.mkdir(parents=True,exist_ok=True)
    base=[str(binary),'pilot',str(family),'-1',str(epsilon),'1.2','unused',str(directory),'-2']
    pilot=subprocess.run(base,check=True,text=True,capture_output=True)
    live,active=map(int,pilot.stdout.strip().split())
    epochs=directory/'epochs.txt'
    epochs.write_text(f'H18 {live+200}\nedge {active-200}\nA21 {active+200}\n')
    (directory/'pilot.json').write_text(json.dumps({'live_sample':live,'active_sample':active})+'\n')
    args=[str(binary),'probe',str(family),'-1',str(epsilon),'1.2',str(epochs),str(directory),str(detailed)]
    subprocess.run(args,check=True)
    # True-bias increments restart from BOOT with a changed admissible root,
    # not from a new root at the word boundary. Each owns its own Live S origin.
    for axis in range(3):
        args[3]=str(axis)
        subprocess.run(args,check=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--eigen',type=Path,default=Path('/usr/include/eigen3'))
    ap.add_argument('--family',type=int,choices=(0,1,2),default=2)
    ap.add_argument('--epsilon',type=float,default=.001)
    ap.add_argument('--detailed-direction',type=int,default=-2)
    ap.add_argument('--skip-build',action='store_true')
    ap.add_argument('--build-only',action='store_true')
    ap.add_argument('--uninstrumented',action='store_true')
    a=ap.parse_args(); a.output=a.output.resolve(); a.output.mkdir(parents=True,exist_ok=True)
    binary=a.output/'shipping-word-probe' if a.skip_build else build(a.output,a.eigen,not a.uninstrumented)
    if not a.build_only:
        run(a.output/f'BIAS{a.family}',binary,a.family,a.epsilon,a.detailed_direction)

if __name__=='__main__': main()
