#!/bin/bash

while getopts "p" OPTION
do
    case $OPTION in
        p)
            print=1
            ;;
    esac
done

g++ random.cpp ex.cpp

./a.out

while getopts "p" OPTION
do
    case $OPTION in
        p)
            jupyter-notebook plot.ipynb
    esac
done

if [[ print -eq 1 ]]
then
    jupyter-notebook plot.ipynb
fi
