#!/bin/bash
# run_tri_lo.sh M aEWM1 tag : LO cross sections of the Dirac-triplet UFO (6 channels) at 13 TeV, nn23lo1
export PATH=/home/claude/mm/env/bin:$PATH
M=$1; AEW=$2; TAG=$3; cd /home/claude/scan/mg
cat > l_tri_$TAG.mg5 <<EOF
launch tri_lo -n tri_$TAG
set nevents 500
set iseed 17
set systematics_program none
set ebeam1 6500
set ebeam2 6500
set MQC $M
set MQ0 $M
set aEWM1 $AEW
EOF
python3 /home/claude/ext/mg5amcnlo/bin/mg5_aMC l_tri_$TAG.mg5 > tri_$TAG.log 2>&1
grep -E "Cross-section" tri_$TAG.log | tail -1
