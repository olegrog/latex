#/bin/bash

d1=2.40014

if ! test -f knfunc.dat; then
    wget http://repository.kulib.kyoto-u.ac.jp/dspace/bitstream/2433/199811/4/knfunc.dat
fi

printf "#$(yes '%8s' | head -13 | xargs)\n" eta Y_0 Y_1/2 Omega_1 -Theta_1 Y_{a4} -Y_{a5} Y_{a6} Omega_3 -Theta_3 Omega_5 -Theta_5 intY_1/2

tail -n +7 knfunc.dat | awk '
NR%5==1 && $1<30 {
    print 1*$1,-$14,-$15,1*$2,-$8,-'$d1'*$15-$16,1*$21,-'$d1'*$15-$20,-$3,1*$9,1*$5,-$11,-$32
}'
