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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;


//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;

  ReadInput >> nspin;

  ReadInput >> J;

  ReadInput >> h;
    
  ReadInput >> metro; // if=1: Metropolis, else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> reload; // if=1: reload configuration from config.final

  ReadInput >> shutup;

  ReadInput >> mag; // magnet field swiched on or off

  if ( shutup == 0 ) {

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;
  cout << "Number of spins = " << nspin << endl;
  cout << "Exchange interaction = " << J << endl;
  cout << "External field = " << h << endl << endl;
  }

  if(metro==1) { if (shutup == 0) cout << "The program performs Metropolis moves" << endl; }
  else { if (shutup == 0) cout << "The program performs Gibbs moves" << endl; }
  if (shutup == 0) cout << "Number of blocks = " << nblk << endl;
  if (shutup == 0) cout << "Number of steps in one block = " << nstep << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

  if(reload == 1) { //Restarting from config.final
    ifstream LoadConf("config.final");

    if ( LoadConf.is_open() ) {
      if (shutup == 0) cout << "Restarting from <<config.final>>." << endl << endl;

      int i = 0;
      int in;
      while (LoadConf >> in) {
        s[i] = in;
        i++;
      }
    }
    else {
      if (shutup == 0) cout << "Cannot find <<config.final>>. Generating random configuration." << endl << endl;
      for (int i=0; i<nspin; ++i)
      {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
      }
    }
    LoadConf.close();
  }
  else {
//initial configuration
  if (shutup == 0) cout << "Generating random configuration." << endl << endl;
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
//  if (shutup == 0) cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      sm = s[o];
      energy_new = Boltzmann(-sm,o);
      energy_old = Boltzmann( sm,o);
      p = exp(-beta*(energy_new-energy_old));

      if(p >= rnd.Rannyu() ) {
        s[o] *= -1;
        accepted += 1;
      }
    }
    else //Gibbs sampling
    {
      energy_up   = 1/( 1 + exp(Boltzmann(-1,o)*2*beta));
      //energy_down = 1/( 1 + exp(Boltzmann(1,o)*2*beta));

      if(energy_up < rnd.Rannyu() ) s[o] = 1;
      else s[o] = -1;
      accepted += 1;
    }
  }
  attempted += nspin;
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, cap = 0.0, m = 0.0, chi = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    if (mag == 0) {
      u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
      cap += pow( -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]) ,2);
    }
    m += s[i];
// INCLUDE YOUR CODE HERE
  }
  walker[iu] = u;
  walker[ic] = cap;
  walker[im] = m;
  walker[ix] = m*m;
// INCLUDE YOUR CODE HERE
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=13, nr = int(temp*8);
    
/*    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
*/
    
    if (mag == 0) {
      Ene.open("output.ene." + to_string(nr).substr(0,4) + ".dat",ios::app);
      stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
      glob_av[iu]  += stima_u;
      glob_av2[iu] += stima_u*stima_u;
      err_u=Error(glob_av[iu],glob_av2[iu],iblk);
      Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Ene.close();
  
      Heat.open("output.hea." + to_string(nr).substr(0,4) + ".dat",ios::app);
      stima_c = beta*beta*( blk_av[ic]/blk_norm/(double)nspin - pow( blk_av[iu]/blk_norm/(double)nspin ,2)); //Heat Capacity
      glob_av[ic]  += stima_c;
      glob_av2[ic] += stima_c*stima_c;
      err_c=Error(glob_av[ic],glob_av2[ic],iblk);
      Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      Heat.close();

      Chi.open("output.chi." + to_string(nr).substr(0,4) + ".dat",ios::app);
      stima_x = blk_av[ix]*beta/blk_norm/(double)nspin; //Magnetic Permeability
      glob_av[ix]  += stima_x;
      glob_av2[ix] += stima_x*stima_x;
      err_x=Error(glob_av[ix],glob_av2[ix],iblk);
      Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.close();
    }

    if (mag == 1) {
      Mag.open("output.mag." + to_string(nr).substr(0,4) + ".dat",ios::app);
      stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
      glob_av[im]  += stima_m;
      glob_av2[im] += stima_m*stima_m;
      err_m=Error(glob_av[im],glob_av2[im],iblk);
      Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.close();
    }

// INCLUDE YOUR CODE HERE

//    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

//  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
