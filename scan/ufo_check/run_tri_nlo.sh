#!/bin/bash
# run_tri_nlo.sh M [aEWM1] : fixed-order NLO QCD cross section (3 of the 6 triplet channels = wino-equivalent total), 13 TeV, nn23nlo
export PATH=/home/claude/mm/env/bin:$PATH
M=$1; AEW=${2:-127.934}; cd /home/claude/scan/mg
cat > l_tri_nlo_$M.mg5 <<EOF
launch tri_nlo -n fo${M}_a${AEW}
fixed_order=ON
set ebeam1 6500
set ebeam2 6500
set MQC $M
set MQ0 $M
set aEWM1 $AEW
set req_acc_FO 0.003
EOF
python3 /home/claude/ext/mg5amcnlo/bin/mg5_aMC l_tri_nlo_$M.mg5 > l_tri_nlo_${M}_a${AEW}.log 2>&1
echo "M=$M aEWM1=$AEW: $(grep 'Total cross section' l_tri_nlo_${M}_a${AEW}.log | tail -1 | sed 's/ \+/ /g') | $(grep -A1 'envelope' l_tri_nlo_${M}_a${AEW}.log | tail -1 | sed 's/ \+/ /g')"
