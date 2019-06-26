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
#include "random2.h"

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

   int N = 1000, M[4] = {1, 2, 10, 100};
   double Sstd[N][4] = {}, Sexp[N][4], Scalo[N][4];
   double lambda = 1, mu = 0, gamma = 1, sigma = 1;
   ofstream out1, out2, out3;
   out1.open("ex2.1.dat");

   for (int j = 0; j < N; j++) {
      for (int i = 0; i < 4; i++) {
        for (int k = 0; k < M[i]; k++) {
            Sstd[j][i] += rnd.Gauss(mu, sigma);
        }
        Sstd[j][i] /= double(M[i]);
      }
      out1 << Sstd[j][0] << " " << Sstd[j][1] << " " << Sstd[j][2] << " " << Sstd[j][3] << endl;
   }
   out1.close();

   out2.open("ex2.2.dat");

   for (int j = 0; j < N; j++) {
      for (int i = 0; i < 4; i++) {
        for (int k = 0; k < M[i]; k++) {
            Sexp[j][i] += rnd.Exp(lambda);
        }
        Sexp[j][i] /= double(M[i]);
      }
      out2 << Sexp[j][0] << " " << Sexp[j][1] << " " << Sexp[j][2] << " " << Sexp[j][3] << endl;
   }
   out2.close();

   out3.open("ex2.3.dat");

   for (int j = 0; j < N; j++) {
      for (int i = 0; i < 4; i++) {
        for (int k = 0; k < M[i]; k++) {
            Scalo[j][i] += rnd.Calo(mu, gamma);
        }
        Scalo[j][i] /= double(M[i]);
      }
      out3 << Scalo[j][0] << " " << Scalo[j][1] << " " << Scalo[j][2] << " " << Scalo[j][3] << endl;
   }
   out3.close();

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
