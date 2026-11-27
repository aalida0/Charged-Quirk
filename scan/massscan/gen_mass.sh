#!/bin/bash
# gen_mass.sh <M_GeV> <nev> : DY (17 -17) and W*-mediated (17 18) samples at mass M, filtered to |eta|<2.5
M=$1; N=${2:-1500}
B=/home/claude/qd_build; S=/home/claude/qd_src
export MAMBA_ROOT_PREFIX=/home/claude/mm/root LD_LIBRARY_PATH=$B/quirkcore
RUN="/home/claude/mm/x/bin/micromamba run -p /home/claude/mm/env"
D=/home/claude/scan/ms/M$M; mkdir -p $D; cd $D
M0=$(python3 -c "print($M-0.3)")
MMIN=$(python3 -c "print($M+$M-0.3+1)")
sed -e '/^HadronLevel:QED/d' -e "s#^PDF:pSet *= *LHAPDF6:.*#PDF:pSet = 20#" $S/gen/pythia_quirk_dy_13TeV.cmnd > dy.cmnd
sed -e "s#^18:m0 = .*#18:m0 = $M0#" -e "s#^34:mMin = .*#34:mMin = $MMIN#" /home/claude/scan/gen/w.cmnd > w.cmnd
( time $RUN $B/quirksim/quirkgen $M $N dy.hepmc3 dy.cmnd 11 35 > dy.log 2>&1 ) 2> dy.time
$RUN python3 $S/gen/prep_gen.py dy.hepmc3 dyf.hepmc3 dy.npz 2.5 > dy.prep 2>&1
( time $RUN $B/quirksim/quirkgen $M $N w.hepmc3 w.cmnd 13 35 > w.log 2>&1 ) 2> w.time
$RUN python3 /home/claude/charged-quirk/scan/prep_gen_pair.py w.hepmc3 wf.hepmc3 w.npz 2.5 > w.prep 2>&1
grep "wrote" dy.log w.log; cat dy.time w.time | grep real; tail -2 dy.prep w.prep
