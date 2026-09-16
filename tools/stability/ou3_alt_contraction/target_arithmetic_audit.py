"""Audit the pinned MCU toolchain without mistaking compilation for a proof.

Run against the verified Arduino 3.3.7 esp-x32/2511 and esp32s3-libs
packages and Arduino Eigen 0.3.2. The report records the actual compiler
options, shipping-header object code and libm identity. It deliberately cannot
qualify target arithmetic, libm accuracy, solver termination or the ALT master.
No download, firmware flash or shipping-source modification is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]
QUALIFICATION = 'OU3_ALT_PINNED_TARGET_ARITHMETIC_AUDIT_V1'
PROBE = r'''
#include "ahrs/Mahony_AHRS.h"
#include "tuner/WavePeriodEstimator.h"
#include "kalman_ou_common/KalmanOUCoreMath.h"
extern "C" void alt_wpe(WavePeriodEstimator* p, float dt, float a) {
    p->update(dt, a);
}
extern "C" void alt_mahony(Mahony_AHRS<float>* p, const float* v, float* angles) {
    p->update(v[0],v[1],v[2],v[3],v[4],v[5],angles,angles+1,angles+2,v[6]);
}
extern "C" void alt_qaxis(float tau, float h, float sigma, float* out) {
    Eigen::Matrix<float,4,4> q;
    ocean_imu::kalman::ou_detail::IntegratedOUChain<float,3>::process_covariance(tau,h,sigma,q);
    Eigen::Map<Eigen::Matrix<float,4,4>> result(out);
    result=q;
}
'''


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(args, **kw):
    result = subprocess.run([str(x) for x in args], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kw)
    if result.returncode:
        raise RuntimeError('target audit command failed: '+shlex.join([str(x) for x in args])+'\n'+result.stderr)
    return result.stdout


def audit(toolchain: Path, sdk: Path, eigen: Path):
    workflow = (ROOT/'.github/workflows/build.yml').read_text()
    for pin in ('esp32:esp32@3.3.7', 'Eigen@0.3.2',
                'esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc',
                '-funroll-loops -fno-finite-math-only'):
        if pin not in workflow:
            raise ValueError('shipping MCU build pin changed: '+pin)
    compiler = toolchain/'bin/xtensa-esp32s3-elf-g++'
    objdump = toolchain/'bin/xtensa-esp32s3-elf-objdump'
    nm = toolchain/'bin/xtensa-esp32s3-elf-nm'
    flags_path = sdk/'flags/cpp_flags'
    flags = shlex.split(flags_path.read_text()) + ['-Os', '-funroll-loops', '-fno-finite-math-only']
    version = run([compiler, '--version']).splitlines()[0]
    if version != 'xtensa-esp-elf-g++ (crosstool-NG esp-14.2.0_20251107) 14.2.0':
        raise ValueError('compiler differs from pinned MCU package')
    if not (eigen/'Eigen/Dense').is_file():
        raise ValueError('Eigen include root must contain Eigen/Dense')
    options = run([compiler, *flags, '-Q', '--help=optimizers'])
    selected = {key: next(line.strip() for line in options.splitlines() if line.lstrip().startswith(key))
                for key in ('-ffp-contract=', '-fassociative-math', '-ffinite-math-only',
                            '-fsigned-zeros', '-frounding-math', '-funsafe-math-optimizations')}
    libm = Path(run([compiler, '-print-file-name=libm.a']).strip())
    if not libm.is_file():
        raise ValueError('target compiler did not resolve libm.a')
    with tempfile.TemporaryDirectory(prefix='ou3-target-audit-') as tmp:
        tmp = Path(tmp)
        source = tmp/'probe.cpp'; source.write_text(PROBE)
        obj = tmp/'probe.o'
        # Standalone numerical-header probe. Use the pinned Eigen sources but
        # do not claim that this includes Arduino's wrapper or final link.
        run([compiler, *flags, '-DEIGEN_NON_ARDUINO', '-DEIGEN_MPL2_ONLY',
             '-ffile-prefix-map='+str(tmp)+'=/ou3-audit',
             '-ffile-prefix-map='+str(ROOT)+'=/ocean-imu',
             '-ffile-prefix-map='+str(toolchain)+'=/toolchain',
             '-ffile-prefix-map='+str(eigen)+'=/eigen',
             '-I'+str(ROOT/'src'), '-I'+str(eigen), '-c', source, '-o', obj])
        disassembly = run([objdump, '-dr', obj])
        undefined = run([nm, '-u', obj])
        instructions = {name: len(re.findall(r'\b'+re.escape(name)+r'\s', disassembly))
                        for name in ('add.s','sub.s','mul.s','madd.s','msub.s','sqrt.s','div.s')}
        object_sha = digest(obj)
    headers = ('ahrs/Mahony_AHRS.h', 'tuner/WavePeriodEstimator.h',
               'kalman_ou_common/KalmanOUCoreMath.h')
    return {
        'qualification': QUALIFICATION,
        'compiler_version': version,
        'compiler_driver_sha256': digest(compiler.resolve()),
        'sdk_cpp_flags_sha256': digest(flags_path),
        'compiler_flags': flags,
        'compiler_options': selected,
        'shipping_headers_sha256': {p:digest(ROOT/'src'/p) for p in headers},
        'eigen_macros_sha256': digest(eigen/'Eigen/src/Core/util/Macros.h'),
        'libm_sha256': digest(libm),
        'libm_multilib_path': str(libm.relative_to(toolchain)),
        'probe_sha256': hashlib.sha256(PROBE.encode()).hexdigest(),
        'probe_object_sha256': object_sha,
        'probe_instruction_counts': instructions,
        'probe_undefined_symbols': sorted(line.split()[-1] for line in undefined.splitlines()),
        'standalone_numerical_header_target_compilation_succeeded': True,
        'probe_only_defines': ['EIGEN_NON_ARDUINO', 'EIGEN_MPL2_ONLY'],
        'complete_Arduino_wrapper_compilation_qualified': False,
        'target_execution_observed': False,
        'whole_firmware_link_resolution_qualified': False,
        'hardware_instruction_rounding_semantics_qualified': False,
        'all_input_libm_error_bounds_qualified': False,
        'all_source_solver_termination_qualified': False,
        'compiler_profile_correspondence_closed': False,
        'storage_search_allowed': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--toolchain', type=Path, required=True)
    parser.add_argument('--sdk', type=Path, required=True)
    parser.add_argument('--eigen', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.toolchain.resolve(), args.sdk.resolve(), args.eigen.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('compiler_version', 'probe_instruction_counts',
                     'standalone_numerical_header_target_compilation_succeeded', 'storage_search_allowed')}))


if __name__ == '__main__':
    main()
