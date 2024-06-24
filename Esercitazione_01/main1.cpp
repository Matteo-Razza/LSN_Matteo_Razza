#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {
    Random rnd;
    initializeRandom(rnd);

    int M = 100000;
    int N = 100; //numero di blocchi
    int L = int(M / N); //numero di estrazioni per ogni blocco
    double ave = 0;
    double ave2 = 0;
    double sum = 0;
    double sum2 = 0;
    double error = 0;
    double accu = 0;
    ofstream flusso_out;

    flusso_out.open("Risultati1.dat");

    if (flusso_out.fail()) {
        cout << endl << "Problema apertura file" << endl;
    }

    for (int i = 0; i < N; i++) {

        for (int j = 0; j < L; j++) {
            accu += rnd.Rannyu();
        }

        ave = accu / L;
        ave2 = pow(ave, 2);

        sum += ave;
        sum2 += ave2;

        if (i == 0) {
            error = 0;
        } else {
            error = sqrt((sum2 / (i + 1) - pow(sum / (i + 1), 2)) / i);
        }

        flusso_out << i + 1 << ' ' << sum / (i + 1) << ' ' << error << endl;

        accu = 0;
    }

    flusso_out.close();

    return 0;
}