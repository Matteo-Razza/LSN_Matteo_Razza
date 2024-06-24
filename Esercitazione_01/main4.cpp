#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {
    Random rnd;
    initializeRandom(rnd);

    int n = 10000;
    double appo41 = 0.;
    double appo42 = 0.;
    double appo43 = 0.;
    double appo51 = 0.;
    double appo52 = 0.;
    double appo53 = 0.;
    double appo61 = 0.;
    double appo62 = 0.;
    double appo63 = 0.;
    ofstream flusso_out4;
    ofstream flusso_out5;
    ofstream flusso_out6;

    flusso_out4.open("Risultati4.dat");
    flusso_out5.open("Risultati5.dat");
    flusso_out6.open("Risultati6.dat");

    if (flusso_out4.fail()) {
        cout << endl << "Problema apertura file" << endl;
    }

    if (flusso_out5.fail()) {
        cout << endl << "Problema apertura file" << endl;
    }

    if (flusso_out6.fail()) {
        cout << endl << "Problema apertura file" << endl;
    }

    for(int j=0; j<n; j++){
        
        for (int i = 0; i < 2; i++){
            appo41 += rnd.Rannyu();
            appo51 += rnd.Exp(1);
            appo61 += rnd.Cauchy(0., 1.);
        }

        for (int i = 0; i < 10; i++){
            appo42 += rnd.Rannyu();
            appo52 += rnd.Exp(1);
            appo62 += rnd.Cauchy(0., 1.);
        }

        for (int i = 0; i < 100; i++){
            appo43 += rnd.Rannyu();
            appo53 += rnd.Exp(1);
            appo63 += rnd.Cauchy(0., 1.);
        }

        flusso_out4 << rnd.Rannyu() << ' ' << appo41/2. << ' ' << appo42/10. << ' ' << appo43/100. << endl;
        flusso_out5 << rnd.Exp(1) << ' ' << appo51/2. << ' ' << appo52/10. << ' ' << appo53/100. << endl;
        flusso_out6 << rnd.Cauchy(0., 1.) << ' ' << appo61/2. << ' ' << appo62/10. << ' ' << appo63/100. << endl;

        appo41 = 0.;
        appo42 = 0.;
        appo43 = 0.;
        appo51 = 0.;
        appo52 = 0.;
        appo53 = 0.;
        appo61 = 0.;
        appo62 = 0.;
        appo63 = 0.;
    }

    flusso_out4.close();
    flusso_out5.close();
    flusso_out6.close();

    return 0;
}