#!/bin/bash -e

proj=$(ls *.geo | head -1 | sed 's/\..*//')

cases=$(ls | grep '^[0-9]' | sort -n)
if test $# -gt 0; then
    cases=$@
fi

time_stat() {
    files="$1/result/*.txt"
    if [ "$(uname)" == "Darwin" ]; then
        stat -l -t '%s' $files | awk '{print $6}' | sort
    else
        stat --format='%Z' $files | sort
    fi
}

humanize() {
    local T=$1
    local D=$((T/60/60/24))
    local H=$((T/60/60%24))
    local M=$((T/60%60))
    local S=$((T%60))
    [[ ! $T ]] && return
    [[ $D > 0 ]] && printf ' %dd' $D
    [[ $H > 0 ]] && printf ' %2dh' $H
    [[ $M > 0 ]] && printf ' %2dm' $M
    [[ $D = 0 && $H = 0 ]] && printf ' %ds' $S
}

fmt=$(printf '%%20s %.0s' {1..4})
date
printf "$fmt\n" 'case' 'percent' 'spend_time' 'remaining_time'
for d in $cases; do
    d=$(basename $d)
    num=$(ls $d/result/m*.txt 2>/dev/null | wc -l)
    if test $num -eq 0; then
        continue
    else
        num=$((num-1))
    fi
    steps=$(grep steps $d/$proj.kep | awk '{ print $2 }' | sed 's/,//')
    macro=$(grep macro $d/$proj.kep | awk '{ print $2 }' | sed 's/,//')
    total=$(echo "($steps-1)/$macro" | bc)
    percent=$(echo "scale=2;100*$num/$total" | bc)
    first_time=$(time_stat $d | head -1)
    last_time=$(time_stat $d | tail -1)
    unset spend_time remaining_time
    if test $first_time -ne $last_time; then
        spend_time=$(($last_time-$first_time))
        if test $total -ne $num; then
            remaining_time=$(echo "scale=0;$spend_time/$num*($total-$num)" | bc)
        fi
    fi
    printf "$fmt\n" $d $percent "$(humanize $spend_time)" "$(humanize $remaining_time)"
done
