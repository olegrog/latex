#!/bin/bash

for f in $(ls log.*); do
    U=${f/*-/}
    cut_x=$(grep cut $f | awk '{ print $3 }' | sed 's/,//')
    cut_r=$(grep cut $f | awk '{ print $6 }' | sed 's/,//')
    ratio=$(grep cut $f | awk '{ print $9 }' | sed 's/,//')
    total=$(grep Total $f | awk '{ print $3 }')
    Ni=$(grep Total $f | sed 's/.*(/(/;s/).*/)/')
    hx=$(grep 'size (x)' $f | awk '{ print $6 }' | sed 's/,//')
    hy=$(grep 'size (y)' $f | awk '{ print $6 }' | sed 's/,//')
    printf "%.1f & %.2f & %.1f & %9s & %d & %.2f & %.2f & %.4f\n" $U $cut_x $cut_r "$Ni" $total $ratio $hx $hy
done
