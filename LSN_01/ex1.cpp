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
#include <iomanip>
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

   rnd.SaveSeed();
   // Here the template ends and the exercice beginns
   // Ex 01.1.1

   int M = 100000, N = 500;
   int L, k, x[N];
   double r[M], mean[M], err[M], sum[M] = {}, s = 0;
   ofstream out1, out2, out3;
   L = M/N;
   for (int i = 0; i < M; i++) {
      r[i] = rnd.Rannyu();
   }

   for (int i = 0; i < N; i++) {
      s = 0;
      for (int j = 0; j < L; j++) {
        k = j+i*L;
        s += r[k];
      }
      sum[i] = s/L;
      x[i] = (i+1)*L;
   }

   for (int i = 0; i < N; i++){
      mean[i] = 0;
      err[i] = 0;
      for (int j = 0; j < i; j++) {
        mean[i] += sum[j];
        err[i] += pow(sum[j],2);
      }
      mean[i] /= (i+1);
      err[i] /= (i+1);
      
      err[i] = sqrt((err[i] - pow(mean[i],2)) / (i+1));
   }
   
   
   out1.open("ex1.1.dat");
   for(int i=0; i<N; i++){
      out1 << x[i] << " " << mean[i] << " " << err[i] << endl;
   }
   out1.close();

   // Ex 01.1.2 (using some if the variables from Ex 1.1)

   double var[N] = {}, errvar[N] = {}, expval = 0.5;
   
   for (int i = 0; i < N; i++) {
      s = 0;
      for (int j = 0; j < L; j++) {
        k = j+i*L;
        s += pow(r[k]-expval,2);
      }
      sum[i] = s/L;
   }

   for (int i = 0; i < N; i++){
      for (int j = 0; j < i; j++) {
        var[i] += sum[j];
        errvar[i] += pow(sum[j],2);
      }
      var[i] /= (i+1);
      errvar[i] /= (i+1);
      
      errvar[i] = sqrt((errvar[i] - pow(var[i],2)) / (i+1));
   }

   out2.open("ex1.2.dat");
   for(int i=0; i<N; i++){
      out2 << x[i] << " " << var[i] << " " << errvar[i] << endl;
   }
   out2.close();

   // Ex 01.1.3

   M = 100;
   N = 1000;
   int J = 100, n = 0;
   double chi2, border[J+1], ra , nM;



   out3.open("ex1.3.dat");
   nM = N/double(M);
   for (int i = 0; i <= M; i++) {
      border[i] = i/double(M);
   }

   for (int p = 0; p < J; p++) {
      s = 0;
      for (int k = 0; k < M; k++) {
        n = 0;
        for (int i = 0; i < N; i++) {
            ra = rnd.Rannyu();
            if ((border[k] <= ra) and (border[k+1] > ra)) {
                n++;
            }
        }
        s += pow((n - nM),2);
      }
      chi2 = s/nM;
      out3 << chi2 << endl;
   }
   out3.close();




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
