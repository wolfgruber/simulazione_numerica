/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <vector>
#include <array>

//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temper,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

// block mean
const int nblock = 100;
std::array <double,2> sum_etot, sum_temp, sum_epot, sum_ekin, sum_pres;
double err_etot, err_temp, err_epot, err_ekin, err_pres;
double block_etot, block_epot, block_ekin, block_temp, block_pres;
std::vector < std::array <double,2> > etot, temp, epot, ekin, pres;
double errsum_etot[nblock], errsum_temp[nblock];
double errsum_epot[nblock], errsum_ekin[nblock];
double errsum_pres[nblock];

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(int);
double Force(int, int);
double Pbc(double);


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
