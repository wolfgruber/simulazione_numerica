/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Ludwig Wolfgruber
//Laboratorio di Simulazone Numerica
//Esercitio 03

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

void Random::launch(void) {
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}


 
int main (int argc, char *argv[]){

    Random rnd;
    rnd.launch();

    // 1) sampling the asset prize directly
    int n = 10000; // number od blocks
    double t = 0, S0 = 100, T = 1, K = 100; // setting the variables to the 
    double r = 0.1, sigma = 0.25, z;  // demandet values
    double C = 0, P = 0, S, sumc = 0, sump = 0;
    double c[n], p[n], errp, errc;
    ofstream out1, out2;  // preparing output

    out1.open("call1.dat");
    out2.open("put1.dat");
    for (int i = 0; i < n; i++) {  // loop over the asset prices
        errc = 0;
        errp = 0;
        z = rnd.Gauss(0,1);
        S = S0 * exp((r-pow(sigma,2)/2)*T+sigma*z*sqrt(T)); // calulating the price in a single step
        c[i] = max(0.,S-K) * exp(-r*T);
        p[i] = max(0.,K-S) * exp(-r*T);
        sumc += c[i];  //summing up for call, put and errors
        sump += p[i];
        C = sumc / (i+1);  // calculating Call(n)
        P = sump / (i+1);  // calculating Put(n)
        for (int k = 0; k < i; k++) { // calculating the error(n)
            errc += pow(C - c[k],2);
            errp += pow(P - p[k],2);
            }
        errc = sqrt(errc/((i+1)*i));
        errp = sqrt(errp/((i+1)*i));

        out1 << C << " " << errc << endl;  //output in files
        out2 << P << " " << errp << endl;
    }
    out1.close();
    out2.close();
    cout << "1) sampling the asset prize directly:" << endl; // output on screen
    cout << "    Call: " << C << " +/- " << errc << endl;
    cout << "    Put:  " << P << " +/- " << errp << endl;

    // 2) sampling discrete path of asset price
    int m = 100; // number of steps for each asset price
    out1.open("call2.dat");
    out2.open("put2.dat");
    sumc = 0;
    sump = 0;
    for (int i = 0; i < n; i++) {
        S = S0; // setting the price to S0 for each loop
        errc = 0, errp = 0;
        for (int j = 0; j < m; j++) { // calculating the price after m teps
            z = rnd.Gauss(0,1);
            S = S * exp((r-pow(sigma,2)/2)*(T/m)+sigma*z*sqrt(T/m));
        }
        c[i] = max(0.,S-K) * exp(-r*T); // calculaiting the call and put profit
        p[i] = max(0.,K-S) * exp(-r*T);
        sumc += c[i]; // summing up and taking the mean of call and put profit
        sump += p[i];
        C = sumc / (i+1);
        P = sump / (i+1);
        for (int k = 0; k < i; k++) { // calculating the error(n)
            errc += pow(C - c[k],2);
            errp += pow(P - p[k],2);
            }
        errc = sqrt(errc/((i+1)*i));
        errp = sqrt(errp/((i+1)*i));
        out1 << C << " " << errc << endl; // output in files
        out2 << P << " " << errp << endl;
    }
    out1.close();
    out2.close();
    cout << "\n2) sampling discrete path of asset price:" << endl; // output on screen
    cout << "    Call: " << C << " +/- " << errc << endl;
    cout << "    Put:  " << P << " +/- " << errp << endl;


    rnd.SaveSeed();
    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
