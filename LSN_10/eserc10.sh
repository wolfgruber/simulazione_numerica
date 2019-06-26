#!/bin/bash

k=2

while getopts "p" OPTION
do
    case $OPTION in
        p)
            print=1
            ;;
    esac
done

while getopts "c:" OPTION
do
	case $OPTION in
		c)
            k=$OPTARG
            ;;
    esac
done

echo "compiling the programm"
g++ random.cpp travel.cpp

sed -i "2s/.*/1/" input.dat
./a.out
cp min.dat cmin.dat
cp conf.dat cconf.dat


sed -i "2s/.*/0/" input.dat
./a.out
cp min.dat smin.dat
cp conf.dat sconf.dat

echo
echo "compiling the paralell program"

mpicxx random.cpp paravel.cpp

mpiexec -np $k a.out

if [[ print -eq 1 ]]
then
    jupyter-notebook plot.ipynb
fi

