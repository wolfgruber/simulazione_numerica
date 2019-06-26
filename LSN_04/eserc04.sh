#!/bin/bash

print=0

while getopts "p:" OPTION
do
    case $OPTION in
        p)
            print=1
            ;;
    esac
done

echo "Compiling the program."
echo
g++ -std=c++11 random.cpp 2MolDyn_NVE.cpp -o esec

cp input.solid input.solid
./solid.sh
cp ave_etot.dat sol_et.dat
cp ave_ekin.dat sol_ek.dat
cp ave_epot.dat sol_ep.dat
cp ave_temp.dat sol_te.dat
cp ave_pres.dat sol_pr.dat

cp input.liquid input.liquid
./liquid.sh
cp ave_etot.dat liq_et.dat
cp ave_ekin.dat liq_ek.dat
cp ave_epot.dat liq_ep.dat
cp ave_temp.dat liq_te.dat
cp ave_pres.dat liq_pr.dat

cp input.gas input.gas
./gas.sh
cp ave_etot.dat gas_et.dat
cp ave_ekin.dat gas_ek.dat
cp ave_epot.dat gas_ep.dat
cp ave_temp.dat gas_te.dat
cp ave_pres.dat gas_pr.dat


if [[ print -eq 1 ]]
then
    jupyter-notebook plot.ipynb
fi

