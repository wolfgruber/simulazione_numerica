/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Ludwig Wolfgruber
// Laboratorio di Simulazione Numerica
// Esercizio 02.1, Montecarlo Integration

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

    int n = 1000;
    double I[n], err[n], sum1 = 0, sum2 = 0, r;
    double pi = M_PI;
    ofstream out1, out2;

    // Integrating with uniform distribution:
    out1.open("int1.dat");
    for (int i = 0; i < n; i++) {
        r = rnd.Rannyu(); // generating uniformly distributed random numbers
        sum1 += pi/2 * cos(pi/2 * r); // summing for the 1st moment
        sum2 += pow(pi/2 * cos(pi/2 * r),2); // summing for the 2nd moment
        I[i] = sum1 / (i+1); // calculating the integral
        err[i] = sum2 / (i+1);
        err[i] = sqrt( (err[i] - pow(I[i],2)) / (i+1) ); // calculating the deviation
        out1 << i+1 << " " << I[i] << " " << err[i] << endl; // output
    }
    out1.close();

    // Integrating with importance sampling ( h(x) = 3/2 (1-x^2) )
    int i = 0;
    double hmax = 3/2, x, c = 0;
    sum1 = 0, sum2 = 0;
    out2.open("int2.dat");
    while (i < n) {
        r = rnd.Rannyu(); // generating two random numbers for the rejection method
        x = rnd.Rannyu();
        if ( 3/2 * (1-pow(x,2)) >= r * hmax) {
        // rejection method to receive random numbers distributed after h(x) = 3/2 (1-x^2)
            sum1 += pi/3 * cos(pi/2 * x)/(1-pow(x,2)); // summing for the 1st moment
            sum2 += pow(pi/3 * cos(pi/2 * x)/(1-pow(x,2)),2); // summing for the 2nd moment
            I[i] = sum1 / (i+1); // calculating the integral
            err[i] = sum2 / (i+1);
            err[i] = sqrt( (err[i] - pow(I[i],2)) / (i+1)); // calculating the deviation
            out2 << i+1 << " " << I[i] << " " << err[i] << " " << x << endl; // output
            i++;
        }
        c++; // counting the total loops
    }
    cout << "Efficiency: " << n/c << endl; // calculating the efficiency of the rejection method

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
