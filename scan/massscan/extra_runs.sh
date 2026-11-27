#!/bin/bash
# sensitivity runs at M=500 GeV: (a) ITk-like card, (b) doublet decay table at the triplet ctau (same Delta)
export MAMBA_ROOT_PREFIX=/home/claude/mm/root LD_LIBRARY_PATH=/home/claude/qd_build/quirkcore
RUN="/home/claude/mm/x/bin/micromamba run -p /home/claude/mm/env"
cd /home/claude/scan/ms
while [ ! -f /home/claude/scan/L3000_c5000_full.finished ]; do sleep 5; done
run() { # dir tag Lambda ctau delta table card hepfile nev
  mkdir -p $1; cat > $1/$2.mac <<EOM
/quirk/mass 500 GeV
/quirk/lambda $3 eV
/quirk/decayCtau $4 mm
/quirk/decayDelta $5 MeV
/quirk/decayTable $6
/qsim/card $7
/qsim/gen/mode hepmc
/qsim/gen/hepmcFile $8
/qsim/gen/quirksOnly true
/qsim/out/printEvery 100000
/run/initialize
/run/beamOn $9
EOM
  (cd $1 && $RUN /home/claude/qd_build/quirksim/quirksim -m $2.mac -t 2 -s 7 -o $2.root > $2.log 2>&1)
}
TAB_T=$PWD/tab_triplet_M500.txt; TAB_D=$PWD/tab_doublet_M500.txt
for P in dy w; do
  for L in 100 300 1000; do
    run out_itk/triplet_${P}_M500 L${L}_c56.2296 $L 56.2296 165.697 $TAB_T /home/claude/qd_src/cards/itk_like.card $PWD/M500/${P}f.hepmc3 1000
    run out_tab/tripletctau_doubtab_${P}_M500 L${L}_c56.2296 $L 56.2296 165.697 $TAB_D /home/claude/qd_src/cards/cms_like.card $PWD/M500/${P}f.hepmc3 1000
  done
done
touch extra.finished
