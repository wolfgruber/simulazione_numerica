
// Ludwig Wolfgruber
// Laboratorio di Simulazione Numerica
// Esercizio 05.1: Variational Monte Carlo

#include <fstream>
#include <cmath>
#include <iostream>
#include "random.h"
#include <vector>

using namespace std;

//some global constants
double pi = 3.14159;
double hbar = 1;
double mass = 1;

double psi(double x, double mu, double sigma) {
    double p = exp(-pow((x-mu),2)/(2*sigma*sigma)) + exp(-pow((x+mu),2)/(2*sigma*sigma));
    return p;
}

double V(double x) {
    return pow(x,4) - 2.5 * pow(x,2);
}

double G(double x, double mu, double sigma) {
    double T, f, g, fp, gp, fpp;
    T = hbar/(2*mass);
    f = pow((x-mu),2)/(2*sigma*sigma);
    g = pow((x+mu),2)/(2*sigma*sigma);
    fp = (x-mu)/(sigma*sigma);
    gp = (x+mu)/(sigma*sigma);
    fpp = pow(sigma,-2);

    T *= (exp(-f) * (fp*fp - fpp) + exp(-g) * (gp*gp - fpp))/psi(x,mu,sigma);
    return V(x) - T;
}
    


vector <double> Energy(double mu, double sigma) {

    Random rnd;
    rnd.launch();

    ofstream out2("psi.dat");

    int n = 40000, m = 1000, c = 0, k = 0, i = 0; // n = 1000000, m = 10000
    int l = n/m, b = 0;
    double x;
    double xold;
    double d, Esum;
    vector <double> Emean, rblock(2,0);
    double p, pold, mo;

    mu = abs(mu);
    sigma = abs(sigma);

    d = sigma + mu;

    xold = mu;
    pold = pow(psi(xold, mu, sigma),2);
    Esum = 0;


    while ( c < n ) {
        x = xold + (rnd.Rannyu()-0.5) * 2 * d;
        p = pow(psi(x, mu, sigma),2);
        
        mo = min(1., p/pold);

        if ( mo > rnd.Rannyu() ) {
            xold = x;
            pold = p;
            c++;
        }

        Esum += G(xold, mu, sigma);

        k++;

        if ( i % m == 0 && i != 0) {
            Emean.push_back(Esum);
            Esum = 0;
            Emean.at(b) /= double(m);

            k = 0;
            
            rblock.at(0) = 0, rblock.at(1) = 0;
            for ( int j = 0; j <= b; j++ ) {
                rblock.at(0) += Emean.at(j);
                rblock.at(1) += pow(Emean.at(j),2);
            }

            rblock.at(0) /= (b+1);

            rblock.at(1) = sqrt( (rblock.at(1)/(b+1) - pow(rblock.at(0),2)) / (b+1) );
            
            b++;
        }

        out2 << xold << endl;
        i++;
    }

    out2.close();
    return rblock;
}


void plot() {
    ofstream Write("plt.dat");
    double mu = 0.5, sigma = 0.35;
    int num = 30;

    for ( int ic = 0; ic < num; ic++) {
        for ( int jc = 0; jc < num; jc++ ) {
            Write << Energy(mu + ic*0.5/double(num), sigma + jc*0.65/double(num)).at(0) << " ";
        }
        Write << endl;
    }
    Write.close();
} 

void min() {
    ifstream Read("input.dat");
    ofstream Write("min.dat");
    double mu, sigma;
    double ds = 0.01, dm = 0.01;
    double nablas, nablam; 
    int num = 30;
    vector <double> E;

    Read >> mu;
    Read >> sigma;
    Read.close();

    for ( int i = 0; i < 50; i++) {

        E = Energy(mu, sigma);
        nablas = (Energy(mu, sigma+ds).at(0) - Energy(mu, sigma-ds).at(0) )/(2*ds);
        nablam = (Energy(mu+dm, sigma).at(0) - Energy(mu-dm, sigma).at(0) )/(2*dm);

        Write << E.at(0) << " " << E.at(1) << " " << mu << " " << sigma << endl;
        cout << "E: " << E.at(0) << " , mu: " << mu << "-" << nablam << " , sigma: " << sigma << "-" << nablas << endl;

        sigma -= nablas*0.01;
        mu -= nablam*0.01;
    }
    cout << "min Energy: " << E.at(0) << " at mu=" << mu << " and sigma=" << sigma << endl;
    Write.close();
} 

int main() {
    //plot();
    min();
} 
            
       









