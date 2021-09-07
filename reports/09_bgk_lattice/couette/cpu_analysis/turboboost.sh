#!/bin/bash

file=log.frequency

: > "$file"_mean
: > "$file"1
: > "$file"2
: > "$file"3

while true; do
    [[ $(ps -u $(whoami) | grep python) ]] || continue
    cat /proc/cpuinfo | grep MHz | awk '{sum += $4} END {print sum/NR}' >> "$file"_mean
    cat /proc/cpuinfo | grep MHz | sort | tail -1 | awk '{print $4}' >> "$file"1
    cat /proc/cpuinfo | grep MHz | sort | tail -2 | head -1 | awk '{print $4}' >> "$file"2
    cat /proc/cpuinfo | grep MHz | sort | tail -3 | head -1 | awk '{print $4}' >> "$file"3
    mean_value=$(cat "$file"_mean | awk '{sum += $1} END {print sum/NR}')
    value1=$(cat "$file"1 | awk '{sum += $1} END {print sum/NR}')
    value2=$(cat "$file"2 | awk '{sum += $1} END {print sum/NR}')
    value3=$(cat "$file"3 | awk '{sum += $1} END {print sum/NR}')
    printf 'Mean = %.0f, first = %.0f, second = %.0f, third = %.0f\n' $mean_value $value1 $value2 $value3
    sleep 0.2
done
