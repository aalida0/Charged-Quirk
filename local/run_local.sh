#!/bin/bash
# run_local.sh <taup|quirk> <mass GeV> <lambda eV> <ngen> <seed> <outdir>
# Same chain as QuirkPheno gen/run_point.sh (quirkgen -> prep_gen.py -> quirksim) with patched Pythia cards.
# Env: QP_WORK (default ~/qp_local, from setup_env.sh), QP_PDFSET (Pythia internal PDF:pSet, default 20 =
# NNPDF3.1 QCD+LUXQED NNLO as=0.118, closest internal set to the cards' LHAPDF6:NNPDF31_nnlo_as_0118).
set -e
M=$1; MASS=$2; L=$3; N=$4; SEED=$5; OUT=$(mkdir -p "$6" && cd "$6" && pwd)
W=${QP_WORK:-$HOME/qp_local}; S=$W/quirkpheno; B=$W/build
export MAMBA_ROOT_PREFIX=$W/root; RUN="$W/x/bin/micromamba run -p $W/env"
export LD_LIBRARY_PATH=$B/quirkcore:$LD_LIBRARY_PATH
CMND=pythia_quirk_dy_13TeV.cmnd; [ "$M" = taup ] && CMND=pythia_taup1e_13TeV.cmnd
# patched card: HadronLevel:QED is unknown to Pythia 8.312; LHAPDF grids not fetchable here -> internal set
sed -e '/^HadronLevel:QED/d' -e "s#^PDF:pSet *= *LHAPDF6:.*#PDF:pSet = ${QP_PDFSET:-20}#" $S/gen/$CMND > $OUT/card.cmnd
T=${M}_m${MASS}_L${L}_s${SEED}
$RUN $B/quirksim/quirkgen $MASS $N $OUT/$T.gen.hepmc3 $OUT/card.cmnd $SEED 35 > $OUT/$T.gen.log 2>&1
$RUN python3 $S/gen/prep_gen.py $OUT/$T.gen.hepmc3 $OUT/$T.genf.hepmc3 $OUT/$T.npz >> $OUT/$T.gen.log 2>&1
NK=$(grep -c "^E " $OUT/$T.genf.hepmc3)
cat > $OUT/$T.mac <<EOM
/quirk/mass $MASS GeV
/quirk/lambda $L eV
/qsim/card $S/cards/cms_like.card
/qsim/gen/mode hepmc
/qsim/gen/hepmcFile $OUT/$T.genf.hepmc3
/qsim/gen/quirksOnly true
/qsim/out/towers true
/qsim/out/printEvery 50
/run/initialize
/run/beamOn $NK
EOM
$RUN $B/quirksim/quirksim -m $OUT/$T.mac -t ${QP_THREADS:-1} -s $SEED -o $OUT/$T.root > $OUT/$T.sim.log 2>&1
echo "$T: generated $N, kept $NK, wrote $OUT/$T.root"
