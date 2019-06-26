// Ludwig Wolfgruber
// Laboratorio di Simulazione Numerica
// Esercizio 09.1, Traveling Salesman

#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>
#include "random.h"

Random rad;

using namespace std;


class Corso { // a class which represents a route through the cities
    public:
    vector<int> it; // the order with which the cities are visited
    double length; // the length of the route

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
vector< vector<double> > InitializeCittaC(int); // initialize cities on a circle
vector< vector<double> > InitializeCittaS(int); // initialize cities inside a square
vector< Corso > InitializeTrips(int, int); // initialize random trips
bool Check(vector< Corso >, int, int); // check if the trips fulfill the bonds
void L1(vector< Corso > &, vector< vector<double> >, int, int); // measuring the length of the trips with L_1-norm
void L2(vector< Corso > &, vector< vector<double> >, int, int); // measuring the length of the trips with L_2-norm
vector< Corso > Sort(vector< Corso >, vector<double>); // sorting the trips according to their length
int Choose(int); // a random number generator of the function x^-6
void OutFinal(Corso, vector< vector< double > >, int); // print the final trip in x,y-coordinates
double Mean(vector< Corso >, int); // calculating the mean length of the shortest trips
vector< Corso > Procreazione(Corso, Corso, int); // creating a new generation member trip
template <typename T> // a function which helps with the sorting
vector<size_t> sort_indexes(const vector<T> &v);


int main() {
    int c = 30, n, l = 0, elite, g, circle;

    ifstream Read("input.dat"); // reading some parameters
    Read >> c;
    Read >> g;
    Read >> elite;
    Read >> n;
    Read >> circle;
    Read.close();

    ofstream out("min.dat"), out2("mean.dat"); // preparing output
    rad.launch(); // launching the Rannyu() random number generator


    int a, b; //indices for parents
    vector< vector<double> > citta; // vector with citie coordinates
    vector< Corso > trips, young, bebe; // trips 
    
    if(circle == 1) citta = InitializeCittaC(c); // initializing the cities 
    else citta = InitializeCittaS(c);

    trips = InitializeTrips(g, c); // initializing the trips

    L1(trips, citta, g, c); // measuring the length of the trips
    sort(trips.begin(), trips.end()); // sorting the trips by length

    for ( int l = 0; l < n; l++ ) { // loop over n generations

        for (int p = 0; p < elite; p++) { // copying the elite (shortest members of the old gen.) to the new generation
            young.push_back(trips.at(p));
        }


        for ( int p = 0; p < (g-elite)/2; p++ ) { // choosing two members of the old gen., mutating and letting them reproduce
            a = Choose(g);
            b = Choose(g);

            bebe = Procreazione(trips.at(a), trips.at(b), c);

            young.push_back(bebe.at(0));
            young.push_back(bebe.at(1));
            bebe.clear();
        }

        trips = young;
        young.clear();

        
        if (!Check(trips, g, c)) { // checking the trips
            cout << "Trips do not fulfil the bonds. The optimisation is stopped" << endl;
            exit(0); // terminating the program if the trips do not fulfill the bonds
        }
        //else cout << "Tutto bene" << endl;

        L1(trips, citta, g, c); // measuring the length of the trips
        sort(trips.begin(), trips.end()); // sorting the trips by length


        if (l % 50 == 0) { // printing some information in cout and to files
            cout << "Generation " << l << ": Min: " << trips.front().length << ", Max: " << trips.back().length << endl;
            out << (*trips.begin()).length << endl;
            out2 << Mean(trips, g) << endl;
            OutFinal(trips.at(0), citta, c);
        }
    }

    out.close(), out2.close();

    OutFinal(trips.at(0), citta, c);
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

vector< Corso > InitializeTrips(int g, int c) {
    vector< Corso > temp;
    vector<int> count, base;
    Corso tt;
    int rint;

    for ( int i = 0; i < c; i++) count.push_back(i);
    
    for ( int j = 0; j < g; j++) {
        base = count;
        for ( int k = 0; k < c; k++) { // creating random trips by randomly choosing numbers out of (0,1,2, ... ,c-1)
            rint = int(rad.Rannyu()*base.size());
            tt.it.push_back(base[rint]);
            base.erase(base.begin()+rint);
        }
        temp.push_back(tt);
        tt.it.clear();
    }
    return temp;
}


bool Check(vector< Corso > trips, int g, int c) {
    bool temp = true;
    vector<int> count, sorted;

    for ( int i = 0; i < c; i++) count.push_back(i);
    
    for ( int j = 0; j < g; j++) { // sorting each trip and comparing it to (0,1,2, ... ,c-1)
        sort(trips[j].it.begin(), trips[j].it.end()); 
        temp = temp and (trips.at(j).it == count);
    }

    return temp;
}

void L1(vector< Corso >& trips, vector< vector<double> > citta, int g, int c) {
    int idx, idx2;
    double t;

    for ( int i = 0; i < g; i++ ) { // looping over the trips in the vector
        idx = trips.at(i).it.at(0);
        idx2 = trips.at(i).it.at(c-1);
        t = sqrt(pow( (citta[idx][0]-citta[idx2][0]),2 )+pow( (citta[idx][1]-citta[idx2][1]),2));
        for ( int j = 1; j < c; j++ ) { // looping over the cities in the order of the current trip and summing the distance
            idx = trips.at(i).it.at(j);
            idx2 = trips.at(i).it.at(j-1);
            t += sqrt(pow( (citta[idx][0]-citta[idx2][0]),2 )+pow( (citta[idx][1]-citta[idx2][1]),2));
        }
        trips[i].length = t;
    }
}

void L2(vector< Corso >& trips, vector< vector<double> > citta, int g, int c) {
    int idx, idx2;
    double t;

    for ( int i = 0; i < g; i++ ) { // looping over the trips in the vector
        idx = trips.at(i).it.at(0);
        idx2 = trips.at(i).it.at(c-1);
        t = sqrt(pow( (citta[idx][0]-citta[idx2][0]),2 )+pow( (citta[idx][1]-citta[idx2][1]),2));
        for ( int j = 1; j < c; j++ ) { // looping over the cities in the order of the current trip and summing the distance
            idx = trips.at(i).it.at(j);
            idx2 = trips.at(i).it.at(j-1);
            t += pow( (citta[idx][0]-citta[idx2][0]),2 )+pow( (citta[idx][1]-citta[idx2][1]),2);
        }
        trips[i].length = t;
    }
}

int Choose(int g) {
    double r = rad.Rannyu();
    return int( ( pow(r,6) ) * (g-1));
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

double Mean(vector< Corso > trips, int g) {
    double temp = 0;
    for (int i = 0; i < g/2; i++) { // looping over the first half of the trips and calculatinf their mean length
        temp += trips.at(i).length;
    }
    return temp*2/g;
}

vector< Corso > Procreazione(Corso a, Corso b, int c) {
    int r = int(rad.Rannyu()*(c-1)), rest = c-r;
    vector< Corso > temp;
    vector< int > sa, sb, ta, tb;
    vector< long unsigned int > ida, idb;

    a.Mutate(0.5, 0.5, 0.1, 0.1, 0.1);
    b.Mutate(0.5, 0.5, 0.1, 0.1, 0.1);

    temp.push_back(a);
    temp.push_back(b);

    a.it.erase(a.it.begin(), a.it.begin()+r);
    b.it.erase(b.it.begin(), b.it.begin()+r);

    ida = sort_indexes(a.it);
    idb = sort_indexes(b.it);
    sort(a.it.begin(), a.it.end());
    sort(b.it.begin(), b.it.end());

    for (int i = 0; i < rest; i++) {
        temp.at(0).it.at( r+idb.at(i) ) = a.it.at(i);
        temp.at(1).it.at( r+ida.at(i) ) = b.it.at(i);
    }
    return temp;
}


template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}



