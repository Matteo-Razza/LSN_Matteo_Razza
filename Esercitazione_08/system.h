#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

class System {

private:
    int _nblocks;         // Number of blocks for block averaging
    int _nsteps;          // Number of simulation steps in each block
    bool _sa;             // Simulated annealing
    double _delta_int;    // Delta of Metropolis algorithm for energy evaluation
    double _delta_ext;    // Delta of Metropolis algorithm for parameters
    double _x0;           // Starting position
    double _mu;           // First parameter
    double _sigma;        // Second parameter
    double _temperature;  // Temperature parameter
    double _beta;         // Inverse temperature parameter
    int _nattempts_int = 0;   // Number of attempted moves in the inner cycle
    int _naccepted_int = 0;   // Number of accepted moves in the inner cycle
    int _nattempts_ext = 0;   // Number of attempted moves in the outer cycle
    int _naccepted_ext = 0;   // Number of accepted moves in the outer cycle
    Random _rnd;          // Random number generator

    double _old_energy;       // Old energy for boltzmann
    double _block_av;         // Block averages of properties
    double _global_av;        // Global averages of properties
    double _global_av2;       // Squared global averages of properties
    double _average;          // Average values of properties

public: // Function declarations
    int get_nbl();                                                                      // Get the number of blocks
    int get_nsteps();                                                                   // Get the number of steps in each block
    void set_beta(double beta);                                                         // Set 1/T
    void set_temperature(double T);                                                     // Set temperature 
    void initialize(bool sa, int blocks, int steps, double delta_int, double x0, vector<double> parameters); // Initialize system 
    void set_delta_ext(double delta_ext);                                               // Set delta Metropolis parameters                       
    void data_blocking(double x0, vector<double> parameters);                           // Data blocking
    double move(double xold,  vector<double> parameters);                               // Move a position
    vector<double> move_parameters(vector<double> old_parameters, int step);            // Move parameters 
    void hamiltonian(vector<double> parameters, double x);                              // Compute hamiltonian
    double psi_2(vector<double> parameters, double x);                                  // Compute wave function
    double boltzmann(double energy);                                                    // Compute Boltzman weight
    void averages(int blk);                                                             // Compute averages of properties
    double error(double acc, double acc2, int blk);                                     // Compute error
    void block_reset();                                                                 // Reset block averages
    void average_reset();                                                               // Reset average accumulators

};

#endif // __System__
