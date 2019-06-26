#!/bin/bash

while getopts "slg" OPTION
do
	case $OPTION in
		s)
			phase=solid
			nstep=1000
			;;
		l)
			phase=liquid
			nstep=1000
			;;
		g)
			phase=gas
			nstep=4000
			;;
	esac
done

g++ -std=c++11 random2.cpp MolDyn_NVE.cpp

cp input2.$phase input2.dat
sed -i "6s/.*/$nstep/" input2.dat
cp config.fcc config.0
echo "Equilibrating the system"
echo
./a.out > /dev/null
cp config.final config.0
cp velocity.final velocity.0

sed -i "8s/.*/0.5/" input2.dat

for ((c=0;c<9;c++))
do
    ./a.out
    cp config.final config.0
    cp velocity.final velocity.0
done

sed -i "6s/.*/100000/" input2.dat

./clean.sh
./a.out 

cp output.gave.1 mol.ga.$phase
cp output.gerr.1 mol.ge.$phase
cp output_epot.dat mol.ep.$phase
cp output_pres.dat mol.pr.$phase
