#!/bin/bash

./clean.sh

g++ -std=c++11 random.cpp Monte_Carlo_NVT.cpp

cp config.fcc config.0
cp input.solid input.dat
sed -i "6s/.*/1/" input.dat
sed -i "7s/.*/250/" input.dat
./a.out
cp config.final config.0
sed -i "6s/.*/1000/" input.dat
sed -i "7s/.*/1/" input.dat
./clean.sh
echo "Computing 10^5 steps in solid phase"
echo
./a.out > /dev/null
cp output.epot.0 outtes.dat
cp output.pres.0 outtps.dat


cp config.fcc config.0
cp input.liquid input.dat
sed -i "6s/.*/1/" input.dat
sed -i "7s/.*/800/" input.dat
./a.out
cp config.final config.0
sed -i "6s/.*/1000/" input.dat
sed -i "7s/.*/1/" input.dat
./clean.sh
echo "Computing 10^5 steps in liquid phase"
echo
./a.out > /dev/null
cp output.epot.0 outtel.dat
cp output.pres.0 outtpl.dat



cp config.fcc config.0
cp input.gas input.dat
sed -i "6s/.*/1/" input.dat
sed -i "7s/.*/4000/" input.dat
./a.out
cp config.final config.0
sed -i "6s/.*/1000/" input.dat
sed -i "7s/.*/1/" input.dat
./clean.sh
echo "Computing 10^5 steps for gas"
echo
./a.out > /dev/null
cp output.epot.0 outteg.dat
cp output.pres.0 outtpg.dat


