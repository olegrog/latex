#!/bin/bash -e

serv=mvs
proj=$(ls *.geo | head -1 | sed 's/\..*//')
dir=${PWD#$HOME/}
[ $# -eq 0 ] && { echo "Usage: ./$(basename $0) <case>"; exit 1; }
case=$1

mkdir $case
scp $serv:$dir/$case/$proj.{geo,kei,msh} $case/

mkdir $case/result
for f in $(ssh $serv ls -t $dir/$case/result | head -5); do
    scp $serv:$dir/$case/result/$f $case/result
done


