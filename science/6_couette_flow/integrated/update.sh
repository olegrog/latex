#!/bin/bash

if test $# -gt 0; then
    files=$@
else
    files="flow qflow shear qflowy pxx pzz temp press diff"
fi

for f in $files; do
    echo $f
    rm -f $f.{pdf,tex}
    make $f.pdf > /dev/null
done
