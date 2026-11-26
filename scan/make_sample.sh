#!/bin/bash
# Drell-Yan Q+Q- sample, M = 500 GeV, 13 TeV, Pythia internal PDF set 20 (see ../docs/LOCAL_RUN.md), |eta| < 2.5 filter.
# usage: make_sample.sh <outdir> [nevents=4000] [seed=11]   (paths below are those of the sandbox build; adjust)
set -e
OUT=$1; N=${2:-4000}; SEED=${3:-11}
B=${QD_BUILD:-/home/claude/qd_build}; S=${QD_SRC:-/home/claude/qd_src}
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-/home/claude/mm/root} LD_LIBRARY_PATH=$B/quirkcore
RUN="${MM:-/home/claude/mm/x/bin/micromamba} run -p ${ENVP:-/home/claude/mm/env}"
mkdir -p $OUT && cd $OUT
sed -e '/^HadronLevel:QED/d' -e "s#^PDF:pSet *= *LHAPDF6:.*#PDF:pSet = 20#" $S/gen/pythia_quirk_dy_13TeV.cmnd > card.cmnd
$RUN $B/quirksim/quirkgen 500 $N gen.hepmc3 card.cmnd $SEED 35 > gen.log 2>&1   # last argument: sigma_z of the beam spot in mm
$RUN python3 $S/gen/prep_gen.py gen.hepmc3 genf.hepmc3 gen.npz 2.5
