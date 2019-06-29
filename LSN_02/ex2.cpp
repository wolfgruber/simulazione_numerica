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
// Esercizio 02.2, Random Walk

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

    int m = 1000, n = 100, r;
    double x[3][m], cx[3][m], walk[n], error[n], sum1, sum2;
    double phi, theta, a = 1;
    ofstream out1("dat2a.dat"), out2;

    // random walk in a discrete cartesian lattice, stepwidth = a
    out1 << 0 << " " << 0 << " " << 0 << endl;
    for (int i = 0; i < n; i++) { // loop over number of steps n
        sum1 = 0;
        for (int j = 0; j < m; j++) { // loop over repetitions
            if (i == 0) { x[0][j] = 0, x[1][j] = 0, x[2][j] = 0; } 
            // setting the origin of each walk equal to zero
            r = int(round(rnd.Rannyu(0.5, 6.5))); // generating random integers in [1,6]
            if (r < 4) { x[r-1][j] += a; }      // walking
            else if (r > 3) { x[r-4][j] -= a; } // walking
            sum1 += pow(x[0][j],2) + pow(x[1][j],2) + pow(x[2][j],2); // computing the squared distance
            //sum2 += pow(pow(x[0][j],2) + pow(x[1][j],2) + pow(x[2][j],2),2); // computing the distance^4 for the error
            sum2 += pow(  pow(x[0][j],2) + pow(x[1][j],2) + pow(x[2][j],2)  ,2);
        }
        error[i] = sum2 / m;
        walk[i] = sum1 / m;
        error[i] = sqrt( (error[i] - pow(walk[i],2)) / (i+1) ); // calculating the error
        out1 << i+1 << " " << sqrt(walk[i]) << " " << sqrt(error[i]) << endl; // output
    }
    out1.close();

    // random walk in a continous lattice with stepwidth a
    out2.open("dat2b.dat");
    out2 << 0 << " " << 0 << " " << 0 << endl;
    for (int i = 0; i < n; i++) { // loop over number of steps n
        sum1 = 0, sum2 = 0;
        for (int j = 0; j < m; j++) { // loop over independet walks
            if (i == 0) { cx[0][j] = 0, cx[1][j] = 0, cx[2][j] = 0; }
            // setting the origin of each walk equal to zero
            phi = rnd.Rannyu() * 2 * M_PI; // generating two random values for the two angles
            theta = acos(1-2*rnd.Rannyu());
            cx[0][j] += a * sin(theta) * cos(phi); // computing the changes on the coordinates
            cx[1][j] += a * sin(theta) * cos(phi);
            cx[2][j] += a * cos(theta);
            sum1 += pow(cx[0][j],2) + pow(cx[1][j],2) + pow(cx[2][j],2); // computing the squared distance
            sum2 += pow( pow(cx[0][j],2) + pow(cx[1][j],2) + pow(cx[2][j],2) ,2);
            //sum2 += pow(pow(cx[0][j],2) + pow(cx[1][j],2) + pow(cx[2][j],2),2);  // computing the distance^4 for the error
        }
        error[i] = sum2 / m;
        walk[i] = sum1 / m;
        error[i] = sqrt( (error[i] - pow(walk[i],2)) / (i+1) ); // calculating the error
        out2 << i+1 << " " << sqrt(walk[i]) << " " << sqrt(error[i]) << endl; // output
    }
    out2.close();

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
