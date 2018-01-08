#!/bin/bash -e

serv=mvs
home=/nethome/aristov
proj=$(ls *.geo | head -1 | sed 's/\..*//')
dir=${PWD#$HOME/}
[ $# -eq 0 ] && { echo "Usage: ./$(basename $0) <case>"; exit 1; }
case=$1
[ $# -eq 2 ] && dir=${dir/\//$2/}
echo $dir

mkdir -p $case/result

rsync -avz $serv:$home/$dir/$case/$proj.{geo,kei,msh,kep} $case/
rsync -avz --exclude=*.bin $serv:$home/$dir/$case/result/ $case/result/


