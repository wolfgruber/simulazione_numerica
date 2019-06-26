#!/bin/bash

plot=0
metro=0
gibbs=0
solut=0
clean=1

while getopts "mgspr" OPTION
do
    case $OPTION in
        m)
            metro=1
            ;;
        g)
            gibbs=1
            ;;
        s)
            solut=1
            metro=1
            gibbs=1
            ;;
        p)
            plot=1
            ;;
        r) #r...data remains
            clean=0
            ;;
    esac
done
if [[ clean -eq 1 ]]
then
    ./clean.sh
fi

g++ -std=c++11 random.cpp exMonte_Carlo_ISING_1D.cpp

if [[ metro -eq 1 ]]
then
    rm -rf out.m.*
    
    sed -i "9s/.*/0/" input.dat # mute -> off
    sed -i "5s/.*/1/" input.dat # metropolis
    sed -i "4s/.*/0.0/" input.dat # h -> 0
    sed -i "10s/.*/0/" input.dat # mag -> off

    for((i=0;i<14;i++))
    do
        awk -v "var=$i" 'NR==1 {$0=0.5+var*0.125} 1' input.dat > temp.dat # set temperature
        cp temp.dat input.dat
	    sed -i "8s/.*/0/" input.dat # restart -> off
	    sed -i "6s/.*/2/" input.dat # set block number
        ./a.out
	    sed -i "8s/.*/1/" input.dat # enable restart and run again
	    sed -i "6s/.*/5/" input.dat # set block number
	    ./a.out
        awk -v "var=$i" 'BEGIN {print "Temperature = " var*0.125+0.5}'

        #extracting the relevant information
        len=$(wc -l output.chi.4.dat | cut -f 1 -d' ')

        sed -n "$len,${len}p" output.ene.$(($i+4)).dat | awk -v "var=$i" '{print 0.5+var*0.125 " " $3 " " $4}' >> out.m.ene.dat
        sed -n "$len,${len}p" output.hea.$(($i+4)).dat | awk -v "var=$i" '{print 0.5+var*0.125 " " $3 " " $4}' >> out.m.hea.dat
        sed -n "$len,${len}p" output.chi.$(($i+4)).dat | awk -v "var=$i" '{print 0.5+var*0.125 " " $3 " " $4}' >> out.m.chi.dat

        sed -i "9s/.*/1/" input.dat # mute -> on

    done

    echo
    sed -i "9s/.*/0/" input.dat # mute -> off
    sed -i "4s/.*/0.02/" input.dat # h -> 0.02
    sed -i "10s/.*/1/" input.dat # mag -> on

    for((i=0;i<14;i++))
    do
        awk -v "var=$i" 'NR==1 {$0=0.5+var*0.125} 1' input.dat > temp.dat # set temperature
        cp temp.dat input.dat
	    sed -i "8s/.*/0/" input.dat # restart -> off
	    sed -i "6s/.*/2/" input.dat # set block number
        ./a.out
	    sed -i "8s/.*/1/" input.dat # enable restart and run again
	    sed -i "6s/.*/5/" input.dat # set block number
	    ./a.out
        awk -v "var=$i" 'BEGIN {print "Temperature = " var*0.125+0.5}'

        #extracting the relevant information
        len=$(wc -l output.mag.4.dat | cut -f 1 -d' ')

        sed -n "$len,${len}p" output.mag.$(($i+4)).dat | awk -v "var=$i" '{print 0.5+var*0.125 " " $3 " " $4}' >> out.m.mag.dat

        sed -i "9s/.*/1/" input.dat # mute -> on

    done

fi

sed -i "9s/.*/0/" input.dat # mute -> off

if [[ gibbs -eq 1 ]]
then
    rm -rf out.g.*

    sed -i "5s/.*/0/" input.dat # gibbs
    sed -i "4s/.*/0.0/" input.dat # h -> 0
    sed -i "10s/.*/0/" input.dat # mag -> off

    for((i=0;i<14;i++))
    do
        awk -v "var=$i" 'NR==1 {$0=0.5+var*0.125} 1' input.dat > temp.dat # set temperature
        cp temp.dat input.dat
	    sed -i "8s/.*/0/" input.dat # restart -> off
	    sed -i "6s/.*/2/" input.dat # set block number
        ./a.out
	    sed -i "8s/.*/1/" input.dat # enable restart and run again
	    sed -i "6s/.*/5/" input.dat # set block number
	    ./a.out
        awk -v "var=$i" 'BEGIN {print "Temperature = " var*0.125+0.5}'

        #extracting the relevant information
        len=$(wc -l output.chi.4.dat | cut -f 1 -d' ')

        sed -n "$len,${len}p" output.ene.$(($i+4)).dat | awk -v "var=$i" '{print 0.5+var*0.125 " " $3 " " $4}' >> out.g.ene.dat
        sed -n "$len,${len}p" output.hea.$(($i+4)).dat | awk -v "var=$i" '{print 0.5+var*0.125 " " $3 " " $4}' >> out.g.hea.dat
        sed -n "$len,${len}p" output.chi.$(($i+4)).dat | awk -v "var=$i" '{print 0.5+var*0.125 " " $3 " " $4}' >> out.g.chi.dat

        sed -i "9s/.*/1/" input.dat # mute -> on

    done

    echo
    sed -i "9s/.*/0/" input.dat # mute -> off
    sed -i "4s/.*/0.02/" input.dat # h -> 0.02
    sed -i "10s/.*/1/" input.dat # mag -> on

    for((i=0;i<14;i++))
    do
        awk -v "var=$i" 'NR==1 {$0=0.5+var*0.125} 1' input.dat > temp.dat # set temperature
        cp temp.dat input.dat
	    sed -i "8s/.*/0/" input.dat # restart -> off
	    sed -i "6s/.*/2/" input.dat # set block number
        ./a.out
	    sed -i "8s/.*/1/" input.dat # enable restart and run again
	    sed -i "6s/.*/5/" input.dat # set block number
	    ./a.out
        awk -v "var=$i" 'BEGIN {print "Temperature = " var*0.125+0.5}'

        #extracting the relevant information
        len=$(wc -l output.mag.4.dat | cut -f 1 -d' ')

        sed -n "$len,${len}p" output.mag.$(($i+4)).dat | awk -v "var=$i" '{print 0.5+var*0.125 " " $3 " " $4}' >> out.g.mag.dat

        sed -i "9s/.*/1/" input.dat # mute -> on

    done

fi

rm output.*

if [[ plot -eq 1 ]]
then
    jupyter-notebook plot.ipnb
fi

