# EcoSIM standalone Ubuntu binary

This bundle contains the standalone, double-precision Release build of EcoSIM
for Ubuntu 24.04 on x86_64, without MPI. The executable is
`bin/ecosim.f90.x`. Input datasets and case working directories are not included.

## Runtime dependencies

On Ubuntu 24.04, open a terminal in the extracted bundle directory and install
the Ubuntu runtime packages listed in `runtime-packages.txt`:

```bash
sudo apt-get update
mapfile -t packages < runtime-packages.txt
sudo apt-get install --no-install-recommends "${packages[@]}"
ldd -r bin/ecosim.f90.x
```

The loader check should report no missing libraries or undefined symbols.
The executable's required GLIBC symbol versions are at most 2.39; this check
does not establish compatibility with other Linux distributions or older Ubuntu
releases.

Use a configured EcoSIM case and its required input files to run simulations.
Check the case's output and restart paths before running the binary in its
working directory.

## Bundle contents

- `SHA256SUMS`: executable checksum; verify with `sha256sum -c SHA256SUMS`.
- `build-info.txt`: source commit, submodules, compiler and build-host details.
- `runtime-packages.txt`: Ubuntu packages owning the resolved runtime libraries.
- `runtime-libraries.txt`: build-host loader verification output.
- `elf-version-info.txt`: executable symbol-version requirements.
- `licenses/`: EcoSIM license and available dependency notices.

Packaging and loader checks verify architecture and runtime linking. They do
not establish the scientific correctness of a simulation.
