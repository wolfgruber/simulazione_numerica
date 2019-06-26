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
// Esercizio 04: Molecular Dynamics

#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <array>
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

double k_boltzmann = 1;

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
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temper;
  cout << "Temperature of the simulation = " << temper << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;

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
        fs = sqrt(3 * temper / sumv2);   // fs = velocity scale factor 
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

     fs = sqrt(3 * temper / sumv2);   // fs = velocity scale factor 
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
  double v, t, vij, p, size;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pres.open("output_pres.dat",ios::app);

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

     p += (pow(1/dr,12) - 0.5 * pow(1/dr,6)); // Pressure

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart;                    //Potential energy per particle
    stima_kin = t/(double)npart;                    //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart;     //Temperature
    stima_etot = (t+v)/(double)npart;                             //Total energy per particle
    stima_pres = rho * stima_temp + 48/(3*vol) * p; // Pressure

    k = nconf % nblock;

    errsum_etot[k] = stima_etot;
    errsum_ekin[k] = stima_kin;
    errsum_epot[k] = stima_pot;
    errsum_temp[k] = stima_temp;
    errsum_pres[k] = stima_pres;

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

    if ( nconf % nblock == 0 ) {

        Epot.open("ave_epot.dat",ios::app);
        Ekin.open("ave_ekin.dat",ios::app);
        Temp.open("ave_temp.dat",ios::app);
        Etot.open("ave_etot.dat",ios::app);
        Pres.open("ave_pres.dat",ios::app);

        sum_epot[0] = 0, sum_epot[1] = 0;
        sum_ekin[0] = 0, sum_ekin[1] = 0;
        sum_etot[0] = 0, sum_etot[1] = 0;
        sum_temp[0] = 0, sum_temp[1] = 0;
        sum_pres[0] = 0, sum_pres[1] = 0;

        for (int i = 0; i < nblock; i++) {
            sum_epot[0] += errsum_epot[i];
            sum_ekin[0] += errsum_ekin[i];
            sum_etot[0] += errsum_etot[i];
            sum_temp[0] += errsum_temp[i];
            sum_pres[0] += errsum_pres[i];
        }

        sum_etot[0] /= (double)nblock;
        sum_ekin[0] /= (double)nblock;
        sum_epot[0] /= (double)nblock;
        sum_temp[0] /= (double)nblock;
        sum_pres[0] /= (double)nblock;

        for (int i = 0; i < nblock; i++) {
            sum_epot[1] += pow(sum_epot[0] - errsum_epot[i],2);
            sum_ekin[1] += pow(sum_ekin[0] - errsum_ekin[i],2);
            sum_etot[1] += pow(sum_etot[0] - errsum_etot[i],2);
            sum_temp[1] += pow(sum_temp[0] - errsum_temp[i],2);
            sum_pres[1] += pow(sum_pres[0] - errsum_pres[i],2);
        }
        sum_epot[1] = sqrt(sum_epot[1]/(double)(nblock-1));
        sum_ekin[1] = sqrt(sum_ekin[1]/(double)(nblock-1));
        sum_etot[1] = sqrt(sum_etot[1]/(double)(nblock-1));
        sum_temp[1] = sqrt(sum_temp[1]/(double)(nblock-1));
        sum_pres[1] = sqrt(sum_pres[1]/(double)(nblock-1));

        epot.push_back(sum_epot);
        ekin.push_back(sum_ekin);
        etot.push_back(sum_etot);
        temp.push_back(sum_temp);
        pres.push_back(sum_pres);

        block_epot = 0, err_epot = 0;
        block_ekin = 0, err_ekin = 0;
        block_etot = 0, err_etot = 0;
        block_temp = 0, err_temp = 0;
        block_pres = 0, err_pres = 0;

        size = (double)epot.size();

        for ( int l = 0; l < (int)size; l++ ) {
            block_epot += epot.at(l)[0];
            err_epot += pow( epot.at(l)[0] ,2);
            block_ekin += ekin.at(l)[0];
            err_ekin += pow( ekin.at(l)[0] ,2);
            block_etot += etot.at(l)[0];
            err_etot += pow( etot.at(l)[0] ,2);
            block_temp += temp.at(l)[0];
            err_temp += pow( temp.at(l)[0] ,2);
            block_pres += pres.at(l)[0];
            err_pres += pow( pres.at(l)[0] ,2);
        }

        block_epot /= size;       
        block_ekin /= size;  
        block_etot /= size;       
        block_temp /= size; 
        block_pres /= size;

        err_epot = sqrt( (err_epot/size - pow(block_epot,2) ) /size );
        err_ekin = sqrt( (err_ekin/size - pow(block_ekin,2) ) /size );
        err_etot = sqrt( (err_etot/size - pow(block_etot,2) ) /size );
        err_temp = sqrt( (err_temp/size - pow(block_temp,2) ) /size );
        err_pres = sqrt( (err_pres/size - pow(block_pres,2) ) /size );

        Epot << block_epot << " " << err_epot << " " << epot.back()[0] << " " << epot.back()[1] << endl;
        Ekin << block_ekin << " " << err_ekin << " " << ekin.back()[0] << " " << ekin.back()[1] << endl;
        Temp << block_temp << " " << err_temp << " " << temp.back()[0] << " " << temp.back()[1] << endl;
        Etot << block_etot << " " << err_etot << " " << etot.back()[0] << " " << etot.back()[1] << endl;
        Pres << block_pres << " " << err_pres << " " << pres.back()[0] << " " << pres.back()[1] << endl;

        Epot.close();
        Ekin.close();
        Temp.close();
        Etot.close();
        Pres.close();
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
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
