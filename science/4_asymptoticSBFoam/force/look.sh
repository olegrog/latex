#!/bin/bash

rename "s/_interpolatedInner//;s/wall/$1/" wall*
./terms.py $1 .5 && open terms-$1.pdf
