#!/bin/bash
# scan.sh "<Lambdas>" "<ctaus>" nev   -> out/L<L>_c<ct>.root (skips existing)
cd /home/claude/scan/out
for L in $1; do for C in $2; do
  T=L${L}_c${C}; [ -f $T.done ] && continue
  ../point.sh $L $C $3 $T && mv ../$T.root ../$T.mac ../$T.log ../${T}_meta.json . 2>/dev/null && touch $T.done
done; done
