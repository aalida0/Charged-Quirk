#!/bin/bash
# ufo_rec_sim.sh <M> "<Lambda list>" : transport of the UFO-based recoil sample (rec_ufo.hepmc3, PDG remapped to 17/18), triplet decay table, as rec_sim.sh
M=$1
export MAMBA_ROOT_PREFIX=/home/claude/mm/root LD_LIBRARY_PATH=/home/claude/qd_build/quirkcore
RUN="/home/claude/mm/x/bin/micromamba run -p /home/claude/mm/env"
cd /home/claude/scan/ms
read DMEV CT < <(awk -v m="M=$M" '$1==m && $2=="triplet"{split($3,a,"=");split($5,b,"=");printf "%.6f %s\n", a[2]*1000, b[2]}' ms_params.txt)
TAB=$PWD/tab_triplet_M$M.txt; OD=out3/triplet_ufo_M$M; mkdir -p $OD
NEV=$(grep -c "^E " M$M/rec_ufo.hepmc3)
for L in $2; do
  TAG=L${L}_c${CT}; [ -f $OD/$TAG.done ] && continue
  cat > $OD/$TAG.mac <<EOM
/quirk/mass $M GeV
/quirk/lambda $L eV
/quirk/decayCtau $CT mm
/quirk/decayDelta $DMEV MeV
/quirk/decayTable $TAB
/qsim/card /home/claude/qd_src/cards/cms_like.card
/qsim/gen/mode hepmc
/qsim/gen/hepmcFile $PWD/M$M/rec_ufo.hepmc3
/qsim/gen/quirksOnly true
/qsim/out/printEvery 100000
/run/initialize
/run/beamOn $NEV
EOM
  (cd $OD && $RUN /home/claude/qd_build/quirksim/quirksim -m $TAG.mac -t 1 -s 7 -o $TAG.root > $TAG.log 2>&1 && touch $TAG.done)
done
