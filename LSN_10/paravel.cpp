// Ludwig Wolfgruber
// Laboratorio di Simulazione Numerica
// Esercizio 10.2, Parallel Simulated Annealing

#include "mpi.h"
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

    MPI::Init(); // initialize MPI
    int size = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();

    int c = 50, nsteps = 10000, l = 0;
    double boltzmann, beta = 1, increment = 10, border = 0.09;
    ofstream out;
    double getmin[size], globact[size];

    int att = 0, acc = 0;
    int run = 1, minrank;
    vector< vector<double> > citta;
    vector<double> temp(2,0);
    double broadcast[2];
    Corso trip, newtrip;

    rad.launch(rank);

    if ( rank == 0 ) { // node 0 prepares the output and initializes the towns
        out.open("paramin.dat");
        cout << "Travelling salesman problem" << endl << "The path through " << c << " cities is minimized with simulated anneiling." << endl;
        cout << size << " nodes are working on the problem" << endl;
        citta = InitializeCittaS(c);
    }

    for (int sync = 0; sync < c; sync++) { // the towns are distributed to the other nodes
        if ( rank == 0 ) {
            broadcast[0] = citta.at(sync).at(0);
            broadcast[1] = citta.at(sync).at(1);
        }
        MPI_Bcast(broadcast,2,MPI_REAL8,0,MPI::COMM_WORLD);
        if ( rank != 0 ) {
            temp.at(0) = broadcast[0];
            temp.at(1) = broadcast[1];
            citta.push_back(temp);
        }
    }

    trip = InitializeTrip(c); // every node works with its own trip
    L1(trip, citta, c);
    newtrip = trip;
    getmin[rank] = trip.length;

    if (!Check(trip, c)) {
        cout << "Trip does not fulfil the bonds. The optimisation is stopped" << endl;
        exit(0);
    }

    while ( run == 1 ) { // the simulated annealing is the same algorithm as in 'travel.cpp'
        newtrip.Mutate(0.5, 0.5, 0.5, 0.5, 0.5);
        L1(newtrip, citta, c);

        if (!Check(newtrip, c)) {
            cout << "Trip does not fulfil the bonds. The optimisation is stopped" << endl;
            exit(0);
        }

        boltzmann = exp(-beta * (newtrip.length-trip.length) );
    
        if ( rad.Rannyu() < boltzmann ) {
            trip = newtrip;
            acc += 1;
        }
        else newtrip = trip; 

        att += 1;
        l += 1;

        if ( trip.length < getmin[rank] ) getmin[rank] = trip.length;

        if (l % nsteps == 0) {
            globact[rank] = acc/double(att);

            MPI_Gather(&getmin[rank],1,MPI_REAL8,getmin,1,MPI_REAL8,0,MPI::COMM_WORLD); // collecting the minimal lengths from the nodes
            MPI_Gather(&globact[rank],1,MPI_REAL8,globact,1,MPI_REAL8,0,MPI::COMM_WORLD); // collecting the acceptance rates from the nodes

            if ( rank == 0 ) { // node 0 keeps track of the border condition and manages the output
                minrank = distance(getmin, max_element(getmin, getmin+size));
                out << getmin[minrank] << endl;
                cout << "Step " << l << ": global min: " << getmin[minrank] << " in node " << minrank <<
                     ", beta = " << beta << ", acception = " << globact[minrank] << endl;
                if ( globact[minrank] < border ) run = 0;
            }

            MPI_Bcast(&run,1,MPI_INTEGER,0,MPI::COMM_WORLD); // node 0 broadcasts the border condition to the others

            acc = 0;
            att = 0;
            beta += increment;
        }
    }
    
    MPI_Bcast(&minrank,1,MPI_INTEGER,0,MPI::COMM_WORLD); // node 0 tells the other nodes who winns

    for ( int sync = 0; sync < c; sync++) { // the winning node sends its trip to node 0
        if ( rank == minrank ) {
            MPI::COMM_WORLD.Send(&trip.it.at(sync),1,MPI::INTEGER,0,sync);
        }
        else if ( rank == 0 ) {
            MPI::COMM_WORLD.Recv(&trip.it.at(sync),1,MPI::INTEGER,minrank,sync);
        }
    }


    if ( rank == 0 ) { // node 0 prints the shortest trip to the output file
        out.close();
        OutFinal(trip, citta, c);
    }

    MPI::Finalize();
}


vector< vector<double> > InitializeCittaC(int c) {
    vector< vector<double> > temp;
    vector<double> tt;
    double r = 1, pi = 3.141597, phi;
    cout << "initializing towns on a circle" << endl;
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
    cout << "initializing towns in a square" << endl;
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
    ofstream out("paraconf.dat");
    for ( int i = 0; i < c; i++) {
        idx = winner.it.at(i);
        out << citta.at(idx).at(0) << " " << citta.at(idx).at(1) << endl;
    }
    idx = winner.it.at(0);
    out << citta.at(idx).at(0) << " " << citta.at(idx).at(1) << endl;
    out.close();
}




