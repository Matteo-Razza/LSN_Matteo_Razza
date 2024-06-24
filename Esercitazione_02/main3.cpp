#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {
  Random rnd;
  initializeRandom(rnd);

  int ntot = 100;  // number of steps for each random walk
  int M = 10000;  // number of simulation of a random walk in 3D always starting
                  // at the origin
  int N = 100;    // number of blocks
  int L = int(M / N);  // number of throws for each block
  vector<double> bin(6, 0);
  double ave = 0;
  double ave2 = 0;
  double sum = 0;
  double sum2 = 0;
  double error = 0;
  double appo = 0;
  double accu = 0;

  ofstream flusso_out;
  ofstream flusso_out2;

  flusso_out.open("Risultati3.dat");
  flusso_out2.open("Risultati4.dat");

  if (flusso_out.fail()) {
    cout << endl << "Problema apertura file" << endl;
  }

  if (flusso_out2.fail()) {
    cout << endl << "Problema apertura file" << endl;
  }

  for (int n = 1; n <= ntot; n++) {
    sum = 0;
    sum2 = 0;

    for (int k = 0; k < N; k++) {
      accu = 0;
      for (int i = 0; i < L; i++) {
        for (int j = 0; j < n; j++) {
          appo = rnd.Rannyu();         // estraggo un numero tra 0 e 1
          bin[(int)(appo * 6)] += 1.;  // incremento il contatore nel bin
                                       // corrispondente alla direzione estratta
        }

        accu += pow(bin[0] - bin[3], 2) + pow(bin[1] - bin[4], 2) +
                pow(bin[2] - bin[5], 2);
        // calcolo il modulo quadro della distanza
        // dall'origine (senza farne la radice) prime tre
        // positive seconde tre negative

        bin.assign(6,
                   0);  // svuoto il vector che teneva le distanze dall'origine
      }

      ave = sqrt(accu / L);  // faccio qui la radice!
      ave2 = pow(ave, 2);

      sum += ave;
      sum2 += ave2;

      if (k == 0) {
        error = 0;
      } else {
        error = sqrt((sum2 / (k + 1) - pow(sum / (k + 1), 2)) / k);
      }

      flusso_out << n << ' ' << k + 1 << ' ' << sum / (k + 1) << ' ' << error
                 << endl;

      if (k == N - 1) {
        flusso_out2 << n << ' ' << sum / (k + 1) << ' ' << error << endl;
      }
    }
  }

  flusso_out.close();
  flusso_out2.close();

  return 0;
}