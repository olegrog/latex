#!/bin/bash

awk -F, '
/^"/ {
    n = 0
    printf "#"
    for (i = 1; i <= NF; i++) {
        if ($i == "\"arc_length\"") { A[0] = i; n++ }
        if ($i == "\"wallMoment:2\"") { A[1] = i; n++ }
    }
}
{
    for (i = 0; i < n; i++) {
        printf "%13s", $A[i]
    }
    printf "\n"
}
' $1 > ${1/0.csv/.raw}

exit 0
