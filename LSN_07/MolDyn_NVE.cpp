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
// Esercizio 07.2: Molecular Dynamics

#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include "random2.h"

using namespace std;

int main(){ 

  Input();             //Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure(nconf);     //Properties measurement
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput, ReadConf, ReadVelo;
  double ep, ek, pr, et, vir, restart;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  /*seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator */

  Random rnd;  // Initialising the random number generator
  rnd.launch();
  
  ReadInput.open("input2.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;
  bin_size = (box/2.0)/(double)nbin;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;

cout << "Restart: " << restart << endl;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Resart from velocity-file
    if ( restart > 0 ) {
      cout << "Read initial velocity from file velocity.0" << endl << endl;
      ReadVelo.open("velocity.0");
      double sumv2 = 0, fs = 1;
      for (int i=0; i<npart; ++i) {
        ReadVelo >> vx[i] >> vy[i] >> vz[i];
        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
      }
      ReadVelo.close();
      if ( restart < 1 ) {
        cout << "Scale velocity to match temp" << endl << endl;
        sumv2 /= (double)npart;
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
      }
      else {
        cout << "Start with temperature from file velocity.0" << endl << endl;
      }
      for (int i=0; i<npart; ++i){
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;

        xold[i] = x[i] - vx[i] * delta;
        yold[i] = y[i] - vy[i] * delta;
        zold[i] = z[i] - vz[i] * delta;
      }
    }

    else {
//Prepare initial velocities
     cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
     double sumv[3] = {0.0, 0.0, 0.0};
     for (int i=0; i<npart; ++i){
        vx[i] = rnd.Rannyu() - 0.5;
        vy[i] = rnd.Rannyu() - 0.5;
        vz[i] = rnd.Rannyu() - 0.5;

        sumv[0] += vx[i];
        sumv[1] += vy[i];
        sumv[2] += vz[i];
     }
     for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
     double sumv2 = 0.0, fs;
     for (int i=0; i<npart; ++i){
        vx[i] = vx[i] - sumv[0];
        vy[i] = vy[i] - sumv[1];
        vz[i] = vz[i] - sumv[2];

        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
     }
     sumv2 /= (double)npart;

     fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
     for (int i=0; i<npart; ++i){
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;

        xold[i] = x[i] - vx[i] * delta;
        yold[i] = y[i] - vy[i] * delta;
        zold[i] = z[i] - vz[i] * delta;
     }
   }
   rnd.SaveSeed();
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(int nconf){ //Properties measurement
  int bin, k;
  double v, t, vij, p;
  double dx, dy, dz, dr;
  double r, vdir, gdir;
  ofstream Epot, Pres, Gofr, Gave, Gerr;

  //Epot.open("output_epot.dat",ios::app);
  //Pres.open("output_pres.dat",ios::app);
  Gofr.open("output.gofr.1",ios::app);
  Gave.open("output.gave.1",ios::app);
  Gerr.open("output.gerr.1",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
     for (int k=0; k<nbin ; k++) {
       if ((dr < (k+1)*bin_size) and (dr >= k*bin_size)) {
         gwalk[k] += 2;
         break;  
       }
     }


     if(dr < rcut){
       p += (pow(1/dr,12) - 0.5 * pow(1/dr,6));
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    /*
    stima_pot = v/(double)npart; //Potential energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_pres = stima_temp + 48/(3*vol) * p; // Pressure
    */
    k = nconf % nblock;
    /*
    errsum_epot[k] = stima_pot;
    errsum_pres[k] = stima_pres;
    */
    for (int n = 0; n < nbin; n++){
      gblk[n] = gblk[n] + gwalk[n];
      gwalk[n] = 0;
    }
    /*
    Epot << stima_pot  << endl;
    Pres << stima_pres << endl;

    Epot.close();
    Pres.close();
    */
    if ( nconf % nblock == 0 ) {
        /*
        Epot.open("ave_epot.dat",ios::app);
        Pres.open("ave_pres.dat",ios::app);

        sum_epot = 0, err_epot = 0;
        sum_pres = 0, err_pres = 0;

        for (int i = 0; i < nblock; i++) {
            sum_epot += errsum_epot[i];
            sum_pres += errsum_pres[i];
        }

        sum_epot /= (double)nblock;
        sum_pres /= (double)nblock;

        for (int i = 0; i < nblock; i++) {
            err_epot += pow(sum_epot - errsum_epot[i],2);
            err_pres += pow(sum_pres - errsum_pres[i],2);
        }
        err_epot = sqrt(err_epot/(double)(nblock-1));
        err_pres = sqrt(err_pres/(double)(nblock-1));
        */
        Gofr << nconf/100;
        Gave << nconf/100;
        Gerr << nconf/100;

        for (int i = 0; i < nbin; i++) {
            r = (i+1) * bin_size;
            vdir = 4*pi*(pow(r+bin_size,3)-pow(r,3))/3;
            gdir = gblk[i]/nblock/(rho * npart * vdir);
            gblk[i] = 0;
    
            gsum[i] += gdir;
            gsum2[i] += gdir*gdir;
            err_gdir = Error(gsum[i],gsum2[i],nconf/100);
        
            Gofr << " " << gdir;
            Gave << " " << gsum[i]/((double)nconf/100);
            Gerr << " " << err_gdir;
        }

        Gofr << endl;
        Gave << endl;
        Gerr << endl;
        /*
        Epot << sum_epot << " " << err_epot << endl;
        Pres << sum_pres << " " << err_pres << endl;

        Epot.close();
        Pres.close();
        */
        Gofr.close();
        Gave.close();
        Gerr.close();

    }

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf, WriteVelo;

  cout << "\nPrint final configuration to file config.final, final velocities to file velocities.final" << endl << endl;
  WriteConf.open("config.final");
  WriteVelo.open("velocity.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelo << vx[i] << " " << vy[i] << " " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelo.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
