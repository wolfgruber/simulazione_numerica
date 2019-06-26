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

g++ random.cpp ex1.cpp -o ex1
./ex1

if [ $? -eq 0 ]
then
    echo "Ercersise 01.1 successfully executed" 
fi

g++ random.cpp ex2.cpp -o ex2
./ex2

if [ $? -eq 0 ]
then
    echo "Ercersise 01.2 successfully executed" 
fi

if [[ print -eq 1 ]]
then
    jupyter-notebook plot.ipynb
fi

