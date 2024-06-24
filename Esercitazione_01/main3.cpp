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

    int M = 100; //divido [0,1] in M sotto intervalli e implemento il test del chi2
    vector<double> bin(M, 0); 
    int n = 10000; //faccio n throws
    int N = 10000; //ripeto il test del chi2 per N volte
    double appo = 0;
    double chi2 = 0;
    ofstream flusso_out;

    flusso_out.open("Risultati3.dat");

    if (flusso_out.fail()) {
        cout << endl << "Problema apertura file" << endl;
    }

    for(int j=0; j<N; j++){  
        
        for (int i = 0; i < n; i++){
            appo = rnd.Rannyu();
            bin[(int)(appo * M)]+=1; //così incremento il contatore nel bin corrispondente
        }
        
        //calcolo il chi2 che è la sommatoria sui sottointervalli di...
        for(int i=0; i<M; i++){ 
            chi2 += pow((bin[i]-double(n/M)),2)/double(n/M); //gli expected sono obv n/M
        }
        
        flusso_out << j+1 << ' ' << chi2 << endl;
        fill(bin.begin(), bin.end(), 0);

        chi2 = 0;
    }

    flusso_out.close();

    return 0;
}