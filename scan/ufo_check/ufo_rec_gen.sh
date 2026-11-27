#!/bin/bash
# ufo_rec_gen.sh <M> : LHE (remapped UFO events, ufo.cmnd) -> quirkgen (Pythia shower) -> pair pT>180, |eta|<2.5 filter -> jets/MET, same as rec_gen.sh
M=$1; PT=${2:-180}
B=/home/claude/qd_build; CQ=/home/claude/charged-quirk/scan
export MAMBA_ROOT_PREFIX=/home/claude/mm/root LD_LIBRARY_PATH=$B/quirkcore
RUN="/home/claude/mm/x/bin/micromamba run -p /home/claude/mm/env"
cd /home/claude/scan/ms/M$M; N=$(grep -c "^<event>" ufo_$M.lhe)
rm -f fifo_u; mkfifo fifo_u
( PAIRPTMIN=$PT $RUN python3 $CQ/prep_gen_pair.py fifo_u rec_ufo.hepmc3 rec_ufo.npz 2.5 > rec_ufo.prep 2>&1 ) &
$RUN $B/quirksim/quirkgen $M $((N-5)) fifo_u ufo.cmnd 21 35 > rec_ufo.log 2>&1
wait; rm -f fifo_u
$RUN python3 $CQ/jetmet_gen.py rec_ufo.hepmc3 jm3_ufo.npz > jm3_ufo.log 2>&1
grep wrote rec_ufo.log; tail -1 rec_ufo.prep; tail -1 jm3_ufo.log
