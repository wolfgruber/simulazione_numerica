#!/bin/bash

while getopts "slg" OPTION
do
	case $OPTION in
		s)
			phase=solid
			nstep=250
			;;
		l)
			phase=liquid
			nstep=800
			;;
		g)
			phase=gas
			nstep=4000
			;;
	esac
done

g++ -std=c++11 random.cpp Monte_Carlo_NVT.cpp
echo "Program compilated"
echo

cp input.$phase input.dat
cp config.fcc config.0
sed -i "6s/.*/1/" input.dat
sed -i "7s/.*/$nstep/" input.dat
echo "Equilibrating the configuration"
echo
./a.out > /dev/null
cp input.$phase input.dat
cp config.final config.0
./clean.sh
./a.out

cp output.gave.0 met.ga.$phase
cp output.gerr.0 met.ge.$phase
cp output.epot.0 met.ep.$phase
cp output.pres.0 met.pr.$phase
