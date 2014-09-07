#!/bin/bash

for f in $@; do
    rename "s/_interpolatedInner//;s/wall/$f/" wall*
    ./terms-new.py $f .5 && open terms-$f.pdf &
done
