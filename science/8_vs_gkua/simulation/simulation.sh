#!/bin/bash -e

num=4

./param.py $num

proj=$(ls *.geo | sed 's/\..*//')
cases=$(ls | grep '^[0-9]')
dir=$(pwd)

if test $# -gt 0; then
    cases=$@
fi

for c in $cases; do
    cd "$dir/$c"
    [ -f log.kes ] && continue
    [ -f log.kes~ ] && continue
    date
    echo "Simulate for $c..."
    if test -f ../m$c.txt; then
        mkdir result
        cp ../m$c.txt result/m0.txt
    fi
    if test -d result; then
        pre="_$(ls | grep 'result$' | sort -r | head -1 | sed 's/r.*//')"
        mv result "$pre"result
        rm -rf VTK
        mv $proj.kei "$pre"$proj.kei
        ../../../tools/nonuniform_ic.py "$pre"$proj.kei $(ls -t "$pre"result/*.txt | head -1) $proj.kei
    elif test $(find . -name '*result' | head -1); then
        pre="$(ls | grep 'result$' | sort -r | head -1 | sed 's/r.*//')"
        rm -rf VTK
    fi
    time mpirun -np $num ../../../src/kes-mpi $proj.kei > log.kes 2>&1
    #echo "$c is done" | mail -s 'Amazon EC2' o.a.rogozin@gmail.com
done

exit 0

cd $dir
./upload_to_micro.sh 172.31.7.210 $cases

echo "Machine is halt" | mail -s 'Amazon EC2' o.a.rogozin@gmail.com
sleep 10
sudo halt
