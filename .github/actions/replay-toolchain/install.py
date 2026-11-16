#!/usr/bin/env python3
"""Extract the qualified compiler without replacing the runner's system packages."""
import argparse
import hashlib
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
import time
import urllib.request


VERSION = '14.2.0-19'
# Official Debian trixie packages, checked against its Packages SHA256 fields.
# Pin the compiler frontend, driver and C++ headers together. GCC's upstream
# version alone does not identify a distribution's compiler backports.
PACKAGES = {
    'cpp-14-x86-64-linux-gnu': 'ef274b5379f5f97fc71619d39ecc84d039d3e184570d207b909c90afaa5d79e0',
    'g++-14-x86-64-linux-gnu': 'c66b009fedd340520279d1af2f91d0668e9f64f0c08fe62179ca1c75e556be58',
    'gcc-14-x86-64-linux-gnu': 'a17ef039f1ba482051c3efb5c2c24070002e60dd0bd09fd456ec31481a11b725',
    'libgcc-14-dev': 'ca6f2d36d96b19b3eb71405b0b80134d8c89380b02204a2512e5c58ceb090628',
    'libstdc++-14-dev': '4b962fac5f1af0bd8b1b3f97e0f47b9d8e9a79e88802143f00ebfbd78ad87a7b',
}


def install(destination):
    if platform.system() != 'Linux' or platform.machine() != 'x86_64':
        raise RuntimeError('the fixed-draw acceptance toolchain requires Linux x86-64')
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    sdk = destination / 'sdk'
    with tempfile.TemporaryDirectory(prefix='packages-', dir=destination) as temporary:
        for package, expected in PACKAGES.items():
            name = f'{package}_{VERSION}_amd64.deb'
            archive = Path(temporary) / name
            url = 'https://deb.debian.org/debian/pool/main/g/gcc-14/' + name
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(url, timeout=45) as response, archive.open('wb') as stream:
                        shutil.copyfileobj(response, stream)
                    break
                except OSError:
                    if attempt == 2:
                        raise
                    time.sleep(attempt + 1)
            actual = hashlib.sha256(archive.read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError(f'SHA256 mismatch for {name}: {actual}')
            subprocess.run(['dpkg-deb', '--extract', str(archive), str(sdk)], check=True)
            print(f'verified {name} {actual}', flush=True)
    binary_dir = destination / 'bin'
    binary_dir.mkdir(exist_ok=True)
    for alias, driver in (('g++', 'g++'), ('c++', 'g++'), ('gcc', 'gcc'), ('cc', 'gcc')):
        link = binary_dir / alias
        link.unlink(missing_ok=True)
        link.symlink_to(sdk / 'usr/bin' / f'x86_64-linux-gnu-{driver}-14')
    compiler = binary_dir / 'g++'
    version = subprocess.check_output([str(compiler), '--version'], text=True)
    if f'Debian {VERSION}' not in version:
        raise RuntimeError(f'unexpected compiler: {version.splitlines()[0]}')
    frontend = Path(subprocess.check_output(
        [str(compiler), '-print-prog-name=cc1plus'], text=True).strip()).resolve()
    if not frontend.is_relative_to(sdk):
        raise RuntimeError(f'compiler frontend escaped the qualified SDK: {frontend}')
    # Exercise driver relocation, C++ headers and the runner's link/runtime
    # dependencies before an expensive simulation build can begin.
    with tempfile.TemporaryDirectory(prefix='probe-', dir=destination) as temporary:
        source = Path(temporary) / 'probe.cpp'
        executable = Path(temporary) / 'probe'
        source.write_text('#include <iostream>\nint main() { std::cout << __VERSION__; }\n')
        subprocess.run([str(compiler), '-std=c++20', str(source), '-o', str(executable)], check=True)
        result = subprocess.check_output([str(executable)], text=True)
        if result != VERSION.split('-')[0]:
            raise RuntimeError(f'unexpected compiled version: {result}')
    print(version.splitlines()[0], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    install(parser.parse_args().destination)
