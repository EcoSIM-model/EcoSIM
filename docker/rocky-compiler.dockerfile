# Use a specific version for stability (matches Rocky Linux 9 lifecycle)
FROM rockylinux:9

# 1. Install DNF plugins and enable the CodeReady Linux Builder (CRB) repository
RUN dnf -y install dnf-plugins-core && \
    dnf config-manager --set-enabled crb

# 2. Run system updates and install core packages (with conflict replacement enabled)
RUN dnf -y update && dnf -y install --allowerasing \
    gcc \
    gcc-c++ \
    gcc-gfortran \
    zip \
    git \
    cmake \
    vim-enhanced \
    gdb \
    automake \
    autoconf \
    libxml2 \
    curl \
    ca-certificates \
    libcurl-devel \
    libxml2-devel \
    libtool \
    python3 \
    && dnf clean all

# Replicate 'python-is-python3' behavior by creating a symlink
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN gcc --version && g++ --version
