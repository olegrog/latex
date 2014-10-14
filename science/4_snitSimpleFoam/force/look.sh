#!/bin/bash

for f in $@; do
    rename "s/_interpolatedInner//;s/wall/$f/" wall*
    ./terms.py $f .5 && open terms-$f.pdf &
done
