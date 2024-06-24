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
  vector<double> r(3, 0);
  double ave = 0;
  double ave2 = 0;
  double sum = 0;
  double sum2 = 0;
  double error = 0;
  double theta = 0;
  double phi = 0;
  double accu = 0;

  ofstream flusso_out;
  ofstream flusso_out2;

  flusso_out.open("Risultati5.dat");
  flusso_out2.open("Risultati6.dat");

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
          theta = acos(1 - 2 * rnd.Rannyu());// estraggo theta e phi UNIFORMEMENTE
          phi = rnd.Rannyu(0, 2*M_PI);
          r[0]+=sin(theta)*cos(phi);
          r[1]+=sin(theta)*sin(phi);
          r[2]+=cos(theta); // coordinate sferiche
        }

        accu += pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2);
        // calcolo il modulo quadro della distanza
        // dall'origine (senza farne la radice)

        r.assign(3,
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

//le incertezze statistiche sono talmente piccole che nel grafico non si vedono...
//vale la pena fare un altro grafico?

//non ripartire dall'origine sempre

//fai il fit

//occhio al discorso dello Jacobiano: guarda lez 1 per campionare angolo solido

/*
Wewanttogeneraterandomvaluesofqandj,suchthatanequal
number of the vectors they describe fall in equal divisions of solid angle. In other words we want to generate uniformly distributed values of j between 0 and 2p (which is easy) and we want to generate values of q between 0 and p distributed according to the frequency function
€
pθ (θ) = 12 sin(θ)
Fθ (θ) = ∫0θ 12 sinθʹdθʹ = 12 −cosθʹθ0 = 12 (1− cosθ)
• Wehave
• Nowwecaninvertthisequationtogiveusq
−1€ #
θ =cos (1−2r) with r ∈$0,1) uniformlydistributed
*/