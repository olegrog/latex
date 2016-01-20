#!/bin/bash -e

cases=$*
functions="A B B1 B2 B3 B4 T0_2 T0_1 T1_2 T1_1 T2_2 T2_1 TT2 TT12 Q3 Q2 QQ3 QQ22"

for f in $functions; do
    echo "Calculate $f"
    for c in $cases; do
        if ! test -f "$c/$f.txt"; then
            ./solve_fredholm.py $f $c | grep I_6 | awk '{print '$c', $3}'
        fi
    done
done

./exact_gamma.py $cases
