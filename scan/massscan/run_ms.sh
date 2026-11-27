#!/bin/bash
# run_ms.sh "<M list>" "<Lambda list>" nev : per-mass Delta m(M), kappa, table, ctau; DY + W*; doublet + triplet
export MAMBA_ROOT_PREFIX=/home/claude/mm/root LD_LIBRARY_PATH=/home/claude/qd_build/quirkcore
RUN="/home/claude/mm/x/bin/micromamba run -p /home/claude/mm/env"
CQ=/home/claude/charged-quirk
cd /home/claude/scan/ms
for M in $1; do
 for REP in doublet triplet; do
  if [ $REP = doublet ]; then T=0.5; Y=0.5; K=1; else T=1; Y=0; K=2; fi
  DM=$(python3 -c "import sys;sys.path.insert(0,'$CQ/ctau');from ctau import delta_m_M;print(delta_m_M($T,$Y,$M))")
  TAB=$PWD/tab_${REP}_M$M.txt
  CT=$(python3 $CQ/quirkdecay/make_decay_table.py $DM $K $M $TAB | awk '/ctau_mm/{printf "%.4f",$2}')
  echo "M=$M $REP dm=$DM GeV ctau=$CT mm" >> ms_params.txt
  DMEV=$(python3 -c "print($DM*1e3)")
  for PROD in dy w; do
   OD=out/${REP}_${PROD}_M$M; mkdir -p $OD
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
/qsim/gen/hepmcFile $PWD/M$M/${PROD}f.hepmc3
/qsim/gen/quirksOnly true
/qsim/out/printEvery 100000
/run/initialize
/run/beamOn $3
EOM
    (cd $OD && $RUN /home/claude/qd_build/quirksim/quirksim -m $TAG.mac -t 2 -s 7 -o $TAG.root > $TAG.log 2>&1 && touch $TAG.done)
   done
  done
 done
done
