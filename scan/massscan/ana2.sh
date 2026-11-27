#!/bin/bash
export MAMBA_ROOT_PREFIX=/home/claude/mm/root LD_LIBRARY_PATH=/home/claude/qd_build/quirkcore
n=$1; M=${n##*_M}
/home/claude/mm/x/bin/micromamba run -p /home/claude/mm/env python3 -I /home/claude/charged-quirk/scan/analyze_ms.py out/$n csv/$n.csv $M > csv/$n.log 2>&1
