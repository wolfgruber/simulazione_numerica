/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
using namespace std;



 
int main (int argc, char *argv[]){

   Random rnd;
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
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   for(int i=0; i<20; i++){
      rnd.Rannyu();
   }


// Begin of the exercise:

   double a = 0, Bp[2], b, d = 1, l = 1, c = 0;
// A, Bp (= B') and B are points in the 2d-plane, where a and b are the y-components of A and B
   int n = 100000, nx = 0, i = 0;
   double pi[n], err[n] = {}, sum1 = 0, sum2 = 0;
   ofstream out;
   out.open("ex3.dat");
   
   while ( i < n ) {
      a = rnd.Rannyu() * d;
// getting a random value for a in the interval [0,d), the x-component of A doesnt matter
      Bp[0] = rnd.Rannyu() - 0.5;
      Bp[1] = rnd.Rannyu() + a - 0.5;
// setting the point Bp in a square around A
        if ( sqrt( pow(Bp[0],2) + pow(Bp[1]-a,2) ) < 0.5 ) {
// we only continue, if Bp lies in a circle with r = 0.5 around A, so we receive a uniform angular distribution
            b = a + l * (Bp[1]-a) / sqrt( pow(Bp[0],2) + pow((Bp[1]-a),2) );
// calculating b as a point on the Vector Bp-A with distance l to the point A
            if ( (b > d) or (b <= 0) ) {
                nx++;
            }
            else {c++;}
// counting the cases when A-B crosses a line
            pi[i] = 2*l*(i+1) / (d*nx);
            if (nx != 0) {
                sum1 = pi[i];
                sum2 = pow(sum1,2);
                err[i] = sqrt( (sum2/(i+1) - pow(sum1/(i+1),2)) );
            }
            out << i+1 << " " << pi[i] << " " << err[i] << endl;
            i++;
// calculating the error and output of data
        }
   }
         
   cout << "Rejected throws: " << c << "\nEfficiency: " << 1-c/(n+c) << endl;   
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
