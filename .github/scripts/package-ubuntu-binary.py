#!/usr/bin/env python3
"""Package the standalone Ubuntu build; fail on unsupported runtime dependencies."""

import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess


def command(*args):
    return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT)


def required_glibc(version_info):
    versions = re.findall(r"\bGLIBC_(\d+(?:\.\d+)+)\b", version_info)
    if not versions:
        raise RuntimeError("No GLIBC version requirements found in executable")
    return max(tuple(map(int, version.split('.'))) for version in versions)


def owning_package(library):
    # Ubuntu's merged-/usr layout can differ from paths recorded by dpkg.
    resolved = str(Path(library).resolve())
    candidates = [str(library), resolved]
    if resolved.startswith('/usr/lib/'):
        candidates.append(resolved.removeprefix('/usr'))
    for candidate in dict.fromkeys(candidates):
        try:
            record = command('dpkg-query', '-S', candidate).splitlines()[0]
            package = record.split(': ', 1)[0]
            if re.fullmatch(r'[a-z0-9][a-z0-9+.-]*(?::amd64)?', package):
                return package
        except subprocess.CalledProcessError:
            pass
    raise RuntimeError(f"Runtime library has no Ubuntu package owner: {library}")


def main():
    if command('uname', '-m').strip() != 'x86_64':
        raise RuntimeError('This package requires an x86_64 build host')
    root = Path.cwd()
    readme = Path(__file__).resolve().with_name('ubuntu_binary.md')
    if not readme.is_file():
        raise RuntimeError(f'Binary package README not found: {readme}')
    source = root / 'build/Linux-x86_64-double-Release/local/bin/ecosim.f90.x'
    if not source.is_file() or not os.access(source, os.X_OK):
        raise RuntimeError(f'Installed EcoSIM executable not found: {source}')
    package = root / 'dist/ecosim-ubuntu24.04-x86_64'
    (package / 'bin').mkdir(parents=True, exist_ok=False)
    binary = package / 'bin/ecosim.f90.x'
    shutil.copy2(source, binary)
    # The build installs absolute CI paths in RPATH. Scientific libraries are
    # static in the default build; remaining libraries must resolve via Ubuntu.
    command('patchelf', '--remove-rpath', str(binary))
    header = command('readelf', '-h', str(binary))
    if 'ELF64' not in header or 'Advanced Micro Devices X86-64' not in header:
        raise RuntimeError('Executable is not a 64-bit x86 ELF')
    versions = command('readelf', '--version-info', str(binary))
    required = required_glibc(versions)
    if required > (2, 39):
        raise RuntimeError(f'Executable requires glibc {required}, newer than 2.39')
    # No LD_LIBRARY_PATH/build tree fallback is allowed during this check.
    env = os.environ.copy()
    env.pop('LD_LIBRARY_PATH', None)
    libraries = subprocess.check_output(
        ['ldd', '-r', str(binary)], text=True, stderr=subprocess.STDOUT, env=env
    )
    if 'not found' in libraries or 'undefined symbol' in libraries:
        raise RuntimeError(libraries)
    paths = re.findall(r'=>\s+(/\S+)', libraries)
    if not paths:
        raise RuntimeError('No resolved runtime libraries found')
    packages = sorted({owning_package(path) for path in paths})
    (package / 'runtime-packages.txt').write_text('\n'.join(packages) + '\n')
    (package / 'runtime-libraries.txt').write_text(libraries)
    (package / 'elf-version-info.txt').write_text(versions)
    metadata = [
        'Target: Ubuntu 24.04, x86_64, glibc <= 2.39 requirement',
        'Mode: standalone, double precision, Release, no MPI',
        f'Maximum required GLIBC symbol version: {".".join(map(str, required))}',
        'Commit: ' + command('git', 'rev-parse', 'HEAD').strip(),
        command('git', 'submodule', 'status', '--recursive'),
        Path('/etc/os-release').read_text(),
        command('getconf', 'GNU_LIBC_VERSION'),
        command('/usr/bin/gfortran-13', '--version'),
        command('/usr/bin/cmake', '--version'),
        command('dpkg-query', '-W', *packages),
        command('file', str(binary)),
    ]
    (package / 'build-info.txt').write_text('\n'.join(metadata))
    shutil.copy2(readme, package / 'README.md')
    licenses = package / 'licenses'
    licenses.mkdir()
    shutil.copy2(root / 'LICENSE.txt', licenses / 'EcoSIM-LICENSE.txt')
    for dependency in ('zlib', 'hdf5', 'netcdf-c', 'netcdf-fortran'):
        dep_root = root / '3rd-partylibs' / dependency
        for name in ('COPYING', 'COPYRIGHT', 'LICENSE', 'LICENSE.txt', 'README'):
            notice = dep_root / name
            if notice.is_file():
                shutil.copy2(notice, licenses / f'{dependency}-{name}')
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    (package / 'SHA256SUMS').write_text(f'{digest}  bin/ecosim.f90.x\n')
    print(f'Packaged {binary}; GLIBC requirement {required}; runtime packages: {packages}')


if __name__ == '__main__':
    main()
