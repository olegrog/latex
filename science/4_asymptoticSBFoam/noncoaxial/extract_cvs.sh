#!/bin/bash

awk -F, '
/^"/ {
    n = 0
    printf "#"
    for (i = 1; i <= NF; i++) {
        if ($i == "\"Points:0\"") { A[0] = i; n++ }
        if ($i == "\"Points:1\"") { A[1] = i; n++ }
        if ($i == "\"U:0\"") { A[2] = i; n++ }
        if ($i == "\"U:1\"") { A[3] = i; n++ }
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
