/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  //int nconf = 1;
  int percent = 0; 
  System SYS;
  SYS.initialize(); // Initialize the System object according to the content of the input files in the ../INPUT/ directory
  SYS.initialize_properties(); // Initialize data members used for measurement of properties
  SYS.block_reset(0); // Reset block accumulators to zero

  //Sometimes it's convenient to do some steps before measuring to garantee equilibration
  int equilibration_steps = 10000;
  cout << "Equilibration is needed:" << endl;
  for(int j=0; j < equilibration_steps; j++){ //loop over steps in a block
    SYS.step();
    // if(j%10 == 0){ //once every 10 loops
    //  SYS.write_XYZ(nconf); //Write actual configuration in XYZ format
    //  nconf++;
    // }
    if(j/(equilibration_steps/100) == percent){
      std::cout << "\rExecution progression: " << percent+1 << "%";
      std::cout.flush();
      percent++;
    }
  }
  cout << endl << "Equilibration completed. Starting measurments..." << endl;

  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      // if(j%10 == 0){ //once every 10 loops
      //   SYS.write_XYZ(nconf); //Write actual configuration in XYZ format
      //   nconf++;
      // }
    }

    percent = i/(SYS.get_nbl()/100.);
    std::cout << "\rExecution progression: " << percent+1 << "%";
    std::cout.flush();

    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }

  cout << endl << "Simulation completed." << endl;
  
  SYS.finalize();

  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
