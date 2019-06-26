#!/bin/bash

while getopts "p" OPTION
do
    case $OPTION in
        p)
            print=1
            ;;
    esac
done

echo "compiling the programm"
g++ -std=c++11 random.cpp travel.cpp

sed -i "5s/.*/1/" input.dat
./a.out
cp min.dat cmin.dat
cp mean.dat cmean.dat
cp conf.dat cconf.dat


sed -i "5s/.*/0/" input.dat
./a.out
cp min.dat smin.dat
cp mean.dat smean.dat
cp conf.dat sconf.dat

if [[ print -eq 1 ]]
then
    jupyter-notebook plot.ipynb
fi
