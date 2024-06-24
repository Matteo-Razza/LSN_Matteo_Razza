#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {
    Random rnd;
    initializeRandom(rnd);

    int M = 100000; //total throws
    int N = 100; //numero di blocchi
    int L = int(M / N); //numero di estrazioni per ogni blocco
    double ave = 0;
    double ave2 = 0;
    double sum = 0;
    double sum2 = 0;
    double error = 0;
    double Si = 0.;
    double S0 = 100.;
    double sigma = 0.25;
    double r = 0.1;
    double K = 100.;
    double C = 0.;
    double T = 1.;
    ofstream flusso_out;

    flusso_out.open("Risultati1.dat");

    if (flusso_out.fail()) {
        cout << endl << "Problema apertura file" << endl;
    }

    for (int i = 0; i < N; i++) {

        for (int j = 0; j < L; j++) {
            Si = S0*exp((r-sigma*sigma/2.)*T+sigma*rnd.Gauss(0.,1.)*sqrt(T));
            C += exp(-r*T) * std::max(0. , Si-K);
        }

        ave = C / L;
        ave2 = pow(ave, 2);

        sum += ave;
        sum2 += ave2;

        if (i == 0) {
            error = 0;
        } else {
            error = sqrt((sum2 / (i + 1) - pow(sum / (i + 1), 2)) / i);
        }

        flusso_out << i + 1 << ' ' << sum / (i + 1) << ' ' << error << endl;

        C=0.;
    }

    flusso_out.close();

    return 0;
}