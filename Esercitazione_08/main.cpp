#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

    int blocks = 100;
    int steps = 1000;
    double delta_int = 2.0;

    vector<double> parameters = {0.616, 0.803}; // sigma, mu
    double x0=0.0;

    System SYS;

    SYS.initialize(0, blocks, steps, delta_int, x0, parameters); // Initialize the System object

    return 0;
}