#!/bin/bash
# Local (no CVMFS/lxplus) build of QuirkPheno with conda-forge. Usage: setup_env.sh [workdir]
# Builds against conda-forge Geant4 (11.4.3 tested; upstream uses LCG_110 / 11.4.1) and Pythia 8.312.
set -e
W=${1:-$HOME/qp_local}; mkdir -p $W; cd $W
[ -x x/bin/micromamba ] || { curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C . bin/micromamba && mkdir -p x && mv bin x/; }
export MAMBA_ROOT_PREFIX=$W/root
x/bin/micromamba create -y -p $W/env -c conda-forge --override-channels \
  geant4 pythia8 hepmc3 cmake ninja make cxx-compiler numpy scipy matplotlib pyyaml lhapdf uproot
[ -d quirkpheno ] || git clone https://github.com/tvami/QuirkPheno quirkpheno   # Apache-2.0, T. A. Vami
x/bin/micromamba run -p $W/env cmake -S quirkpheno -B build -G Ninja -DCMAKE_PREFIX_PATH=$W/env
x/bin/micromamba run -p $W/env ninja -C build -j${NJOBS:-2}
echo "built in $W/build"
