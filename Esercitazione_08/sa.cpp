#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

    ofstream coutp("parameters.dat");

    double percent = 0;

    int blocks = 100;
    int steps = 10000;
    double temperature = 1; 
    int external_steps = 1000;

    double delta_int = 2.0; // sampling x
    double delta_ext = 1.; // for the parameters
    vector<double> parameters = {0.6, 0.8}; // sigma, mu
    double x0 = 0.0;

    System SYS;

    SYS.initialize(true, blocks, steps, delta_int, x0, parameters); // Initialize the System object

    SYS.average_reset();
    SYS.block_reset(); // Reset block accumulators to zero

    SYS.set_delta_ext(delta_ext);

    for(int k=0; k < external_steps; k++){

        SYS.block_reset();
        SYS.average_reset();

        SYS.set_beta(1./temperature+1.5*k);

        parameters = SYS.move_parameters(parameters, k); // Perform a MC step

        coutp <<  parameters[0] << " " << parameters[1] << endl;

        // percent = double(k)/double(external_steps)*100;
        // std::cout << "\rExecution progression: " << percent << "%";
        // std::cout.flush();
    }

    cout << endl;

    coutp.close();

    return 0;
}