#!/bin/bash -e

proj=$(ls *.geo | head -1 | sed 's/\..*//')

cases="$(ls | grep '^[0-9]' | sort -n) $(ls | grep '.*_[0-9]' | sort -n)"
if test $# -gt 0; then
    cases=$@
fi

for d in $cases; do
    if test ! -f $d/log.kes && test ! -d $d/result; then
        continue
    fi
    (
        cd $d
        if test ! -f log.emptyFoam; then
            echo " --- Handling $d..."
            ../Allclean
            ../Allrun
        fi
    )
done

#make -j8 && open pics/*.pdf
