#!/bin/bash
# rec_gen.sh <M_GeV> <dy|w> <nev> [pairptmin=180] : recoil-enriched sample: generate <nev> events through a FIFO, keep (|eta|<2.5 as before) only
# events with generator-level pair pT > pairptmin, then jets/MET of the kept events, and (DY) the LO coupling re-weighting.
M=$1; P=$2; N=$3; PT=${4:-180}
B=/home/claude/qd_build; S=/home/claude/qd_src; CQ=/home/claude/charged-quirk/scan
export MAMBA_ROOT_PREFIX=/home/claude/mm/root LD_LIBRARY_PATH=$B/quirkcore
RUN="/home/claude/mm/x/bin/micromamba run -p /home/claude/mm/env"
D=/home/claude/scan/ms/M$M; cd $D
[ $P = dy ] && { PROC=11; CMND=dy.cmnd; } || { PROC=13; CMND=w.cmnd; }
rm -f fifo_$P; mkfifo fifo_$P
( PAIRPTMIN=$PT $RUN python3 $CQ/prep_gen_pair.py fifo_$P rec_$P.hepmc3 rec_$P.npz 2.5 > rec_$P.prep 2>&1 ) &
$RUN $B/quirksim/quirkgen $M $N fifo_$P $CMND $((PROC+10)) 35 > rec_$P.log 2>&1
wait; rm -f fifo_$P
python3 $CQ/jetmet_gen.py rec_$P.hepmc3 jm3_$P.npz > jm3_$P.log 2>&1
[ $P = dy ] && $RUN python3 $CQ/reweight_dy.py rec_dy.hepmc3 rec_dy_rw.npz $M > rec_dy_rw.log 2>&1
grep wrote rec_$P.log; cat rec_$P.prep | tail -1; tail -1 jm3_$P.log
