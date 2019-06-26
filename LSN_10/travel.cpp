// Ludwig Wolfgruber
// Laboratorio di Simulazione Numerica
// Esercizio 10.1, Simulated Annealing

#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>
#include "random.h"

Random rad;

using namespace std;


class Corso {
    public:
    vector<int> it;
    double length;

    bool operator<(Corso b) { return length < b.length; }

    //mutations:
    void Shift();
    void PairPerm();
    void MShift();
    void MPerm();
    void Invers();
    void Mutate(double, double, double, double, double);
};

void Corso::Shift() {
    int r = int(rad.Rannyu()*it.size());
    rotate(it.begin(), it.begin()+r, it.end());
}

void Corso::PairPerm() {
    int t, r = int(rad.Rannyu()*(it.size()-1));
    t = it.at(r);
    it.at(r) = it.at(r+1);
    it.at(r+1) = t;
}

void Corso::MShift() {
    int r = int(rad.Rannyu()*it.size());
    int m = int(rad.Rannyu()*(it.size()-r));
    int i = int(rad.Rannyu()*(it.size()-m));
    vector<int> t;

    t.assign(it.begin()+r, it.begin()+r+m);
    it.erase(it.begin()+r, it.begin()+r+m);
    it.insert(it.begin()+i, t.begin(), t.end());
}

void Corso::MPerm() {
    int m = int(rad.Rannyu()*(it.size()/2-1));
    int r1 = int(rad.Rannyu()*(it.size()-2*m));
    int r2 = int(rad.Rannyu()*(it.size()-2*m-r1)+m+r1);
    vector<int> t1, t2;

    t1.assign(it.begin()+r1, it.begin()+r1+m);
    t2.assign(it.begin()+r2, it.begin()+r2+m);

    it.erase(it.begin()+r2, it.begin()+r2+m);
    it.erase(it.begin()+r1, it.begin()+r1+m);

    it.insert(it.begin()+r1, t2.begin(), t2.end());
    it.insert(it.begin()+r2, t1.begin(), t1.end());
}

void Corso::Invers() {
    int r = int(rad.Rannyu()*it.size());
    int m = int(rad.Rannyu()*(it.size()-r));
    vector<int> t;

    t.assign(it.begin()+r, it.begin()+r+m);

    for ( int i = 0; i < m; i++) {
        it.at(r+i) = t.at(m-i-1);
    }
}

void Corso::Mutate(double p1, double p2, double p3, double p4, double p5) {

    if (rad.Rannyu() > p1) {
        Invers();
    }

    if (rad.Rannyu() > p2) {
        Shift();
    }

    if (rad.Rannyu() > p3) {
        MShift();
    }

    if (rad.Rannyu() > p4) {
        MPerm();
    }

    if (rad.Rannyu() > p5) {
        PairPerm();
    }

}


//functions:
vector< vector<double> > InitializeCittaC(int);
vector< vector<double> > InitializeCittaS(int);
Corso InitializeTrip(int);
bool Check(Corso, int);
void L1(Corso &, vector< vector<double> >, int);
void L2(Corso &, vector< vector<double> >, int);
void OutFinal(Corso, vector< vector< double > >, int);



int main() {
    int c, nsteps, l = 0, circle;
    double boltzmann, beta, increment, border;

    ifstream Read("input.dat"); // reading parameters from the input file
    Read >> c; // number of cities
    Read >> circle;
    Read >> beta; // initial temperature, T = 1/beta
    Read >> increment; //beta-increment per round
    Read >> border; // lower acceptance limit
    Read >> nsteps; // number of steps per round
    Read.close();

    ofstream out("min.dat"); // output file
    rad.launch();

    int att = 0, acc = 0; // attempt/accept counters
    bool run = true;
    vector< vector<double> > citta;
    Corso trip, newtrip;
    
    if(circle == 1) citta = InitializeCittaC(c); // initializing the towns
    else citta = InitializeCittaS(c);

    trip = InitializeTrip(c); // initialize the trip
    L1(trip, citta, c); // measuring the length
    newtrip = trip; 

    if (!Check(trip, c)) { // checking if the trip fulfills the bonds
        cout << "Trip does not fulfil the bonds. The optimisation is stopped" << endl;
        exit(0);
    }

    while ( run ) { // loop runs as long as acceptance < border
        newtrip.Mutate(0.5, 0.5, 0.5, 0.5, 0.5); // mutate newtrip randomly
        L1(newtrip, citta, c); // measure it

        if (!Check(newtrip, c)) { // check the mutated trip
            cout << "Trip does not fulfil the bonds. The optimisation is stopped" << endl;
            exit(0);
        }
        //else cout << "Tutto bene" << endl;

        boltzmann = exp(-beta * (newtrip.length-trip.length) ); // calculate the probability of a change to newtrip
    
        if ( rad.Rannyu() < boltzmann ) {
            trip = newtrip; // execute change
            acc += 1;
        }
        else newtrip = trip; // reset newtrip

        att += 1;
        l += 1;

        if (l % nsteps == 0) { // check the border condition and print the min. trip every nsteps
            cout << "Step " << l << ": Min: " << trip.length << ", beta = " << beta << ", acception = " << acc/double(att) << endl;
            out << trip.length << endl;
            if ( acc/double(att) < border ) run = false;
            acc = 0;
            att = 0;
            beta += increment;
        }
    }

    out.close();

    OutFinal(trip, citta, c);
}


vector< vector<double> > InitializeCittaC(int c) {
    vector< vector<double> > temp;
    vector<double> tt;
    double r = 1, pi = 3.141597, phi;
    cout << endl << "initializing towns on a circle" << endl;
    for ( int i = 0; i < c; i++) {
        phi = rad.Rannyu()*2*pi;
        tt.push_back(cos(phi)*r);
        tt.push_back(sin(phi)*r);
        temp.push_back(tt);
        tt.clear();
    }
    return temp;
}

vector< vector<double> > InitializeCittaS(int c) {
    vector< vector<double> > temp;
    vector<double> tt;
    int a = 1;
    cout << endl << "initializing towns in a square" << endl;
        for ( int i = 0; i < c; i++) {
        tt.push_back(rad.Rannyu()*a);
        tt.push_back(rad.Rannyu()*a);
        temp.push_back(tt);
        tt.clear();
    }
    return temp;
}

Corso InitializeTrip(int c) {
    vector<int> base;
    Corso tt;
    int rint;

    for ( int i = 0; i < c; i++) base.push_back(i);

    for ( int k = 0; k < c; k++) {
        rint = int(rad.Rannyu()*base.size());
        tt.it.push_back(base.at(rint));
        base.erase(base.begin()+rint);
    }

    return tt;
}


bool Check(Corso trip, int c) {
    vector<int> count, sorted;
    bool temp;

    for ( int i = 0; i < c; i++) count.push_back(i);
    
    sort(trip.it.begin(), trip.it.end());
    temp = (trip.it == count);

    return temp;
}

void L1( Corso & trips, vector< vector<double> > citta, int c) {
    int idx, idx2;
    double t;

    idx = trips.it.at(0);
    idx2 = trips.it.at(c-1);
    t = sqrt(pow( (citta[idx][0]-citta[idx2][0]),2 )+pow( (citta[idx][1]-citta[idx2][1]),2));
    for ( int j = 1; j < c; j++ ) {
        idx = trips.it.at(j);
        idx2 = trips.it.at(j-1);
        t += sqrt(pow( (citta[idx][0]-citta[idx2][0]),2 )+pow( (citta[idx][1]-citta[idx2][1]),2));
    }

    trips.length = t;
}

void L2( Corso & trips, vector< vector<double> > citta, int c) {
    int idx, idx2;
    double t;


    idx = trips.it.at(0);
    idx2 = trips.it.at(c-1);
    t = sqrt(pow( (citta[idx][0]-citta[idx2][0]),2 )+pow( (citta[idx][1]-citta[idx2][1]),2));
    for ( int j = 1; j < c; j++ ) {
        idx = trips.it.at(j);
        idx2 = trips.it.at(j-1);
        t += pow( (citta[idx][0]-citta[idx2][0]),2 )+pow( (citta[idx][1]-citta[idx2][1]),2);
    }
    trips.length = t;
}

    
void OutFinal(Corso winner, vector< vector< double > > citta, int c) {
    int idx;
    ofstream out("conf.dat");
    for ( int i = 0; i < c; i++) {
        idx = winner.it.at(i);
        out << citta.at(idx).at(0) << " " << citta.at(idx).at(1) << endl;
    }
    idx = winner.it.at(0);
    out << citta.at(idx).at(0) << " " << citta.at(idx).at(1) << endl;
    out.close();
}




