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
    double x = 0; //posizione del punto medio dell'ago rispetto alla riga di sx
    double theta = 0; //angolo acuto tra la verticale e l'ago
    double l=1; //lunghezza dell'ago
    double d=2; //distanza tra le righe verticali

    double Nhit = 0;
    
    double prob = 0;
    double mypi = 3;
    double mypi2 = 0;

    double sum = 0;
    double sum2 = 0;
    double error = 0;
    
    ofstream flusso_out;

    flusso_out.open("Risultati7.dat");

    if (flusso_out.fail()) {
        cout << endl << "Problema apertura file" << endl;
    }

    for (int i = 0; i < N; i++) {

        for (int j = 0; j < L; j++) {
            x = rnd.Rannyu(0, d);
            theta = rnd.Rannyu(0, mypi/2); //"if possible do not use pi to evaluate pi"
                                           //forse è intelligente guessare mypi=3 e sostituire mypi con quello determinato alla fine dell'algoritmo  
                                           //man mano sarà più preciso...

            if(x<=0.5*l*sin(theta) || (2-x)<=0.5*l*sin(theta)){  
                Nhit += 1;
            }
        }

        prob = Nhit/L;
        mypi = (2*l)/(d*prob);
        mypi2 = pow(mypi, 2);

        sum += mypi;
        sum2 += mypi2;

        if (i == 0) {
            error = 0;
        } else {
            error = sqrt((sum2 / (i + 1) - pow(sum / (i + 1), 2)) / i);
        }

        flusso_out << i + 1 << ' ' << sum / (i + 1) << ' ' << error << endl;

        Nhit = 0;
    }

    flusso_out.close();

    return 0;
}