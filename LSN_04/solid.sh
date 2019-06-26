#!/bin/bash

./clean.sh

cp input.solid input.dat
sed -i "6s/.*/2000/" input.dat
#sed -i "7s/.*/500/" input.dat

cp config.fcc config.0

./esec

cp config.final config.0
cp velocity.final velocity.0

sed -i "8s/.*/0.5/" input.dat

for ((c=0;c<9;c++))
do
    ./esec
    cp config.final config.0
    cp velocity.final velocity.0
done

./clean.sh
sed -i "6s/.*/30000/" input.dat
#sed -i "7s/.*/500/" input.dat

./esec


