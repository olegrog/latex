#!/bin/bash -e

proj=$(ls *.geo | head -1 | sed 's/\..*//')

cases=$(ls | grep '^[0-9]' | sort -n)
if test $# -gt 0; then
    cases=$@
fi

date
printf "%7s %7s\n" 'work' 'percent'
for d in $cases; do
    d=$(basename $d)
    num=$(ls $d/result/*.txt 2>/dev/null | wc -l)
    if ! test $num; then
        continue
    fi
    steps=$(grep steps $d/$proj.kep | awk '{ print $2 }' | sed 's/,//')
    macro=$(grep macro $d/$proj.kep | awk '{ print $2 }' | sed 's/,//')
    total=$(echo $steps/$macro | bc)
    percent=$(echo "scale=2;100*$num/$total" | bc)
    printf "%7s %.1f\n" $d $percent
done
