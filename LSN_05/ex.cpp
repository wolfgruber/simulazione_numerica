
// Ludwig Wolfgruber
// Laboratorio di Simulazione Numerica
// Esercizio 05

#include <fstream>
#include <cmath>
#include <iostream>
#include "random.h"

using namespace std;

double psi1(double x, double y, double z) {
    double psi, r, phi, theta, pi = M_PI, a0 = 1;
    r = sqrt( x*x + y*y + z*z );

    psi = sqrt( pow(a0,3) / pi ) * exp(-a0*r);
    return pow(psi,2);
}

double psi2(double x, double y, double z) {
    double psi, r, phi, theta, pi = M_PI, a0 = 1;
    r = sqrt( x*x + y*y + z*z );
    theta = acos(z/r);

    psi = sqrt(2*pow(a0,5)/pi)/8 * r * exp(-a0*r/2) * cos(theta);
    return pow(psi,2);
}

int main() {

    Random rnd;
    rnd.launch();

    ofstream out1("s-orbit.dat"), out2("p-orbit.dat");
    ofstream outr1("s-rad.dat"), outr2("p-rad.dat");

    int n = 10000, m = 100, c = 0, k = 0; // n = 1000000, m = 10000
    int l = n/m, b = 0;
    double x, y, z;
    double xold, yold, zold;
    double d = 2.45, rsum[m], rerr[l], rmean[l];
    double rblock, errblock;
    double p, pold, mo;

    xold = 0, yold = 0, zold = 0;
    pold = psi1(xold,yold,zold);

    for ( int i = 0; i < n; i++) {
        x = xold + (rnd.Rannyu()-0.5) * d;
        y = yold + (rnd.Rannyu()-0.5) * d;
        z = zold + (rnd.Rannyu()-0.5) * d;
        p = psi1(x,y,z);
        
        mo = min(1., p/pold);

        if ( mo > rnd.Rannyu() ) {
            xold = x, yold = y, zold = z;
            pold = p;
            c++;
        }

        rsum[k] = sqrt( pow(xold,2) + pow(yold,2) + pow(zold,2) );
        k++;

        if ( i % m == 0 && i != 0) {
            rmean[b] = 0, rerr[b] = 0;
            for ( int j = 0; j < m; j++ ) {
                rmean[b] += rsum[j];
            }
            rmean[b] /= double(m);
            for ( int j = 0; j < m; j++ ) {
                rerr[b] += pow(rsum[j] - rmean[b], 2);
            }
            rerr[b] = sqrt( rerr[b] / double(m-1) );
            k = 0;
            
            rblock = 0, errblock = 0;
            for ( int j = 0; j <= b; j++ ) {
                rblock += rmean[j];
                errblock += pow(rmean[j],2);
            }

            rblock /= (b+1);

            errblock = sqrt( (errblock/(b+1) - pow(rblock,2)) / (b+1) );
            
            outr1 << rblock << " " << errblock << endl;
            b++;
        }

        out1 << xold << " " << yold << " " << zold << endl;
    }

    out1.close(), outr1.close();
    cout << "Acceptance at s-orbital: " << c/double(n) << ", with d = " << d << endl;

    c = 0, b = 0, k = 0;

    xold = 0, yold = 0, zold = 5;
    pold = psi1(xold,yold,zold);

    d = 5.96;

    for ( int i = 0; i < n; i++) {
        x = xold + (rnd.Rannyu()-0.5) * d;
        y = yold + (rnd.Rannyu()-0.5) * d;
        z = zold + (rnd.Rannyu()-0.5) * d;
        p = psi2(x,y,z);
        
        mo = min(1., p/pold);

        if ( mo > rnd.Rannyu() ) {
            xold = x, yold = y, zold = z;
            pold = p;
            c++;
        }

        rsum[k] = sqrt( pow(xold,2) + pow(yold,2) + pow(zold,2) );
        k++;

        if ( i % m == 0 && i != 0) {
            rmean[b] = 0, rerr[b] = 0;
            for ( int j = 0; j < m; j++ ) {
                rmean[b] += rsum[j];
            }
            rmean[b] /= double(m);
            for ( int j = 0; j < m; j++ ) {
                rerr[b] += pow(rsum[j] - rmean[b], 2);
            }
            rerr[b] = sqrt( rerr[b] / double(m-1) );
            k = 0;
            
            rblock = 0, errblock = 0;
            for ( int j = 0; j <= b; j++ ) {
                rblock += rmean[j];
                errblock += pow(rmean[j],2);
            }

            rblock /= (b+1);

            errblock = sqrt( (errblock/(b+1) - pow(rblock,2)) / (b+1) );
            
            outr2 << rblock << " " << errblock << endl;
            b++;
        }

        out2 << xold << " " << yold << " " << zold << endl;
    }

    cout << "Acceptance at p-orbital: " << c/double(n) << ", with d = " << d << endl;

    out2.close(), outr2.close();
}















