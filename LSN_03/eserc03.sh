#!/bin/bash

print=0

while getopts "p" OPTION
do
    case $OPTION in
        p)
            print=1
            ;;
    esac
done

g++ random.cpp ex.cpp -o ex
./ex

if [ $? -eq 0 ]
then
    echo
    echo "Ercersise successfully executed" 
fi


if [[ print -eq 1 ]]
then
    jupyter-notebook plot.ipynb
fi

