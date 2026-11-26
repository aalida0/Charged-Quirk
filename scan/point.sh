#!/bin/bash
# point.sh <Lambda_eV> <ctau_mm> <nev> <tag> [mass=500] [table]
L=$1; CT=$2; N=$3; TAG=$4; M=${5:-500}; TAB=${6:-/home/claude/charged-quirk/quirkdecay/macros/triplet_dm167.txt}
export MAMBA_ROOT_PREFIX=/home/claude/mm/root LD_LIBRARY_PATH=/home/claude/qd_build/quirkcore
cd /home/claude/scan
cat > $TAG.mac <<EOM
/quirk/mass $M GeV
/quirk/lambda $L eV
/quirk/decayCtau $CT mm
/quirk/decayDelta ${DELTA:-167} MeV
/quirk/decayTable $TAB
/qsim/card /home/claude/qd_src/cards/cms_like.card
/qsim/gen/mode hepmc
/qsim/gen/hepmcFile /home/claude/scan/genf.hepmc3
/qsim/gen/quirksOnly true
/qsim/out/printEvery 100000
/run/initialize
/run/beamOn $N
EOM
/home/claude/mm/x/bin/micromamba run -p /home/claude/mm/env /home/claude/qd_build/quirksim/quirksim -m $TAG.mac -t 2 -s 7 -o $TAG.root > $TAG.log 2>&1
