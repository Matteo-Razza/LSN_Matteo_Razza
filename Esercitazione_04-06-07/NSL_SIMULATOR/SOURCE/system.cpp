/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"

using namespace std;
using namespace arma;

///////////////////////////////////////////////////////////////////

int System :: get_nbl(){ // Get number of blocks
  return _nblocks;
}

int System :: get_nsteps(){ // Get number of steps
  return _nsteps;
}

///////////////////////////////////////////////////////////////////

void System :: initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#    N_BLOCK:      ACCEPTANCE:" << endl;
  couta.close();

  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){
      input >> _sim_type;
      if(_sim_type > 1){
        input >> _J;
        input >> _H;
      }
      if(_sim_type > 3){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      if(_sim_type == 0)      coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"  << endl;
      else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"         << endl;
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl;
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl;
    } else if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "TEMP" ){
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << endl;
    } else if( property == "NPART" ){
      input >> _npart;
      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){ 
        _particle(i).initialize();
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho;
      _side.resize(_ndim);
      _halfside.resize(_ndim);
      double side = pow(_volume, 1.0/3.0);
      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside=0.5*_side; //armadillo*scalar
      coutf << "SIDE= ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i];
      }
      coutf << endl;
    } else if( property == "R_CUT" ){
      input >> _r_cut;
      coutf << "R_CUT= " << _r_cut << endl;
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration();
  this->initialize_velocities();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void System :: read_configuration(){ // Read configuration from a .xyz file in directory ../INPUT/CONFIG/

  string file_input;

  if(_sim_type == 0){
    file_input = "config.xyz";
  } else if(_sim_type == 2){
    file_input = "config.ising";  
  } else if(_sim_type == 3){
    file_input = "config.ising";  
  } else if(_sim_type == 1){
    file_input = "config.xyz";
  }

  ifstream cinf;
  cinf.open("../INPUT/CONFIG/" + file_input);

  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != _npart){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & " + file_input + " not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    for(int i=0; i<_npart; i++){
      cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0)); //Scale coordinates from file into box dimensions
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else cerr << "PROBLEM: Unable to open INPUT file " + file_input<< endl;
  cinf.close();
  if(_restart and _sim_type > 1){
    int spin;
    cinf.open("../INPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++){
      cinf >> spin;
      _particle(i).setspin(spin);
    }
    cinf.close();
  }
  return;
}

void System :: initialize_velocities(){ // Initialize velocities
  if(_restart and _sim_type==0){ //If restart read previous velocities
    ifstream cinf;
    cinf.open("../INPUT/CONFIG/velocities.in");
    if(cinf.is_open()){
      double vx, vy, vz;
      for(int i=0; i<_npart; i++){
        cinf >> vx >> vy >> vz;
        _particle(i).setvelocity(0,vx);
        _particle(i).setvelocity(1,vy);
        _particle(i).setvelocity(2,vz);
      }
    } else cerr << "PROBLEM: Unable to open INPUT file velocities.in"<< endl;
    cinf.close();
  } else { // Only if _sim_type == 0???
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim);
    sumv.zeros();
    for (int i=0; i<_npart; i++){
      vx(i) = _rnd.Gauss(0.,sqrt(_temp)); //Maxwell-Boltzmann distributions
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      sumv(0) += vx(i); //Compute drift velocities
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }
    for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart);
    double sumv2 = 0.0, scalef;
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0); //Subtract drift velocity per particle
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i); //But we have to restore original temperature
    }
    sumv2 /= double(_npart);
    scalef = sqrt(3.0 * _temp / sumv2);   //Velocity scale factor 
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef); //Scale velocities
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
  }
  if(_sim_type == 0){ //_xold initialization for Verlet algorithm 
  double xold, yold, zold;
    for (int i=0; i<_npart; i++){
      // xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0); //True=current position, False=previous position
      // yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      // zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);

      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta + 0.5*_delta*_delta*force(i,0), 0); // Higher order expansion
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta + 0.5*_delta*_delta*force(i,1), 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta + 0.5*_delta*_delta*force(i,2), 2);

      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }
  return;
}

void System :: initialize_properties(){ // Initialize data members used for measurement of properties

  string property;
  int index_property = 0;
  _nprop = 0;

  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;
  _measure_temp     = false;

  ifstream input("../INPUT/properties.dat");
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){
        ofstream coutp("../OUTPUT/potential_energy.dat");
        coutp << "#      BLOCK:      ACTUAL_PE:           PE_AVE:             ERROR:" << endl;
        coutp.close();
        _nprop++;
        _index_penergy = index_property;
        _measure_penergy = true;
        index_property++;
        _vtail = 8.*M_PI*_rho*(1./(9.*pow(_r_cut, 9))-1./(3.*pow(_r_cut,3))); // TO BE FIXED IN EXERCISE 7
        //cout << "_vtail: " << _vtail << endl;
      } else if( property == "KINETIC_ENERGY" ){
        ofstream coutk("../OUTPUT/kinetic_energy.dat");
        coutk << "#      BLOCK:      ACTUAL_KE:           KE_AVE:             ERROR:" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){
        ofstream coutt("../OUTPUT/total_energy.dat");
        coutt << "#      BLOCK:      ACTUAL_TE:           TE_AVE:             ERROR:" << endl;
        coutt.close();
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        _vtail = 8.*M_PI*_rho*(1./(9.*pow(_r_cut, 9))-1./(3.*pow(_r_cut,3)));
        // Ha senso anche qui nel caso in cui properties contempli TOTAL_ENERGY e non POTENTIAL_ENERGY
        index_property++;
      } else if( property == "TEMPERATURE" ){
        ofstream coutte("../OUTPUT/temperature.dat");
        coutte << "#      BLOCK:      ACTUAL_T:            T_AVE:             ERROR:" << endl;
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      } else if( property == "PRESSURE" ){
        ofstream coutpr("../OUTPUT/pressure.dat");
        coutpr << "#      BLOCK:      ACTUAL_P:            P_AVE:             ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        _ptail = 32.*M_PI*_rho*(1./(9.*pow(_r_cut, 9))-1./(6.*pow(_r_cut, 3))); // TO BE FIXED IN EXERCISE 7
        //cout << "_ptail: " << _ptail << endl;
      } else if( property == "GOFR" ){
        ofstream coutgr("../OUTPUT/gofr.dat");
        coutgr << "#    DISTANCE:     AVE_GOFR:            ERROR:" << endl;
        coutgr.close();
        input>>_n_bins;
        _nprop+=_n_bins;
        _bin_size = (_halfside.min() )/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property+= _n_bins;
      } else if( property == "MAGNETIZATION" ){
        ofstream coutpr("../OUTPUT/magnetization.dat");
        coutpr << "#      BLOCK:      ACTUAL_M:            M_AVE:             ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){
        ofstream coutpr("../OUTPUT/specific_heat.dat");
        coutpr << "#      BLOCK:      ACTUAL_CV:           CV_AVE:             ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){
        ofstream coutpr("../OUTPUT/susceptibility.dat");
        coutpr << "#      BLOCK:      ACTUAL_X:            X_AVE:             ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;

  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}

///////////////////////////////////////////////////////////////////

void System :: step(){ // Perform a simulation step
  if(_sim_type == 0) this->Verlet();  // Perform a MD step
  else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly choosen particle
  _nattempts += _npart; //update number of attempts performed on the system
  return;
}

void System :: Verlet(){ // Verlet algorithm scheme
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //Force acting on particle i
    _fx(i) = this->force(i,0);
    _fy(i) = this->force(i,1);
    _fz(i) = this->force(i,2);
  }
  for(int i=0; i<_npart; i++){ //Verlet integration scheme
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0); //Unitary mass in Lennard-Jones units
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);

    // _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    // _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    // _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));

    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,true), 0)/_delta + 0.5*_delta*_fx(i));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,true), 1)/_delta + 0.5*_delta*_fy(i));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,true), 2)/_delta + 0.5*_delta*_fz(i));

    _particle(i).acceptmove(); // xold = xnew
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart;
  return;
}

double System :: force(int i, int dim){ // Force on particle i
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt( dot(distance,distance) ); //inner product between armadillo
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  return f;
}

void System :: move(int i){ // Propose a MC move for particle i
  if(_sim_type == 3){         //Gibbs sampler for Ising
    // TO BE FIXED IN EXERCISE 6
    double acceptance;
    acceptance = 1./(1.+exp(-2.0*_beta*(_J*(_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin())+_H)));
    if(_rnd.Rannyu() < acceptance ){
      _particle(i).setspin(1);
    }else{
      _particle(i).setspin(-1);
    }
    _naccepted++;
  } else {                    // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}

bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() * 
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}

double System :: Boltzmann(int i, bool xnew){ // Boltzmann weight
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}

double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
  return position - _side(i) * rint(position / _side(i));
}

int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

///////////////////////////////////////////////////////////////////

void System :: measure(){ // Measure properties
  _measurement.zeros();
  int bin;
  vec distance;
  distance.resize(_ndim);
  double dr=0.0;
  double penergy_temp=0.0; // temporary accumulator for potential energy
  double kenergy_temp=0.0; // temporary accumulator for kinetic energy
  double tenergy_temp=0.0; // temporary accumulator for total energy
  double temp_temp=0.0;    // temporary accumulator for temperature
  double magn_temp=0.0;    // temporary accumulator for magnetization
  double susc_temp=0.0;    // temporary accumulator for susceptibility
  double virial=0.0;
 
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  if (_measure_penergy or _measure_tenergy or _measure_gofr or _measure_pressure or _measure_cv) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) );
        
        // GOFR ... TO BE FIXED IN EXERCISE 7
        if(dr<_halfside.min() and _measure_gofr) _measurement(_index_gofr+int(dr/_bin_size))+=2;
        
        if(dr < _r_cut){
          if(_measure_penergy or _measure_tenergy or _measure_cv) penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6); // POTENTIAL ENERGY
          if(_measure_pressure) virial += 1.0/pow(dr,12) - 0.5/pow(dr,6); // VIRIAL 
        }
      }
    }
  }
  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy or _measure_tenergy or _measure_cv){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    if (_measure_penergy) _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy or _measure_tenergy or _measure_temp or _measure_pressure or _measure_cv){
    for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
    kenergy_temp /= double(_npart);
    if (_measure_kenergy) _measurement(_index_kenergy) = kenergy_temp;
  }
  // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
  if (_measure_tenergy or _measure_cv){
    if (_sim_type < 2){
      tenergy_temp = kenergy_temp + penergy_temp;
      if (_measure_tenergy) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    } else { 
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      tenergy_temp /= double(_npart);
      if (_measure_tenergy) _measurement(_index_tenergy) = tenergy_temp;
    }
  }  
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp or _measure_pressure){
    temp_temp = (2.0/3.0) * kenergy_temp;
    if (_measure_temp) _measurement(_index_temp) = temp_temp;
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    _measurement(_index_pressure) = _ptail + _rho*temp_temp + 1.0/(3.0*_volume)*48.0*virial;
  }
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if (_measure_magnet){
    double s_i;
    for (int i=0; i<_npart; i++){
      s_i = double(_particle(i).getspin());
      magn_temp += s_i;
    }
    _measurement(_index_magnet) = magn_temp;
  }
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if (_measure_cv){
    _measurement(_index_cv) = pow(tenergy_temp*double(_npart),2); // Only for H=0
  }
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if (_measure_chi){
    double s_i;
    for (int i=0; i<_npart; i++){
      s_i = double(_particle(i).getspin());
      susc_temp += s_i;
    }
    _measurement(_index_chi) = _beta*pow(susc_temp,2); // Only for H=0
  }

  _block_av += _measurement; //Update block accumulators

  return;
}

void System :: write_XYZ(int nconf){ // Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << setprecision(8) << _particle(i).getposition(0,true)          // x
            << setw(16) << setprecision(8) << _particle(i).getposition(1,true)          // y
            << setw(16) << setprecision(8) << _particle(i).getposition(2,true) << endl; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}

void System :: averages(int blk){ // Data blocking

  ofstream coutf;
  double average, sum_average, sum_ave2;

  if(_measure_gofr){
    for(int r=0; r<_n_bins; r++){
      _block_av(_index_gofr + r) /= _rho * double(_npart)* 4.0/3.0 * M_PI * ( pow((r+1)*_bin_size,3) - pow(r*_bin_size,3) ); //normalization
    }
  }

  _average     = _block_av / double(_nsteps);

  if (_measure_cv){
    _average(_index_cv) = pow(_beta, 2)*(_average(_index_cv)- pow(_npart*_average(_index_tenergy), 2));
  }

  _global_av  += _average;
  _global_av2 += _average % _average; // % -> element-wise multiplication

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << setw(10) << blk
          << setw(20) << setprecision(10) << average
          << setw(20) << setprecision(10) << sum_average/double(blk)
          << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << setw(10) << blk
          << setw(20) << setprecision(10) << average
          << setw(20) << setprecision(10) << sum_average/double(blk)
          << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TOTAL ENERGY //////////////////////////////////////////////////////////////
  if (_measure_tenergy){
    coutf.open("../OUTPUT/total_energy.dat",ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    coutf << setw(10) << blk
          << setw(20) << setprecision(10) << average
          << setw(20) << setprecision(10) << sum_average/double(blk)
          << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp){
    coutf.open("../OUTPUT/temperature.dat",ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << setw(10) << blk
          << setw(20) << setprecision(10) << average
          << setw(20) << setprecision(10) << sum_average/double(blk)
          << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    coutf.open("../OUTPUT/pressure.dat",ios::app);
    average  = _average(_index_pressure);
    sum_average = _global_av(_index_pressure);
    sum_ave2 = _global_av2(_index_pressure);
    coutf << setw(10) << blk
          << setw(20) << setprecision(10) << average
          << setw(20) << setprecision(10) << sum_average/double(blk)
          << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // GOFR //////////////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 7
  if (_measure_gofr and blk==_nblocks){
    coutf.open("../OUTPUT/gofr.dat",ios::app);
    for(int r=0; r<_n_bins; r++){
      average  = _average(_index_gofr+r);
      sum_average = _global_av(_index_gofr+r);
      sum_ave2 = _global_av2(_index_gofr+r);
      coutf << setw(10) << r*_bin_size
            << setw(20) << setprecision(10) << sum_average/double(blk)
            << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    }
    coutf.close();
  }
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if (_measure_magnet){
    coutf.open("../OUTPUT/magnetization.dat",ios::app);
    average  = _average(_index_magnet);
    sum_average = _global_av(_index_magnet);
    sum_ave2 = _global_av2(_index_magnet);
    coutf << setw(10) << blk
          << setw(20) << setprecision(10) << average
          << setw(20) << setprecision(10) << sum_average/double(blk)
          << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if (_measure_cv){
    coutf.open("../OUTPUT/specific_heat.dat",ios::app);
    average  = _average(_index_cv);
    sum_average = _global_av(_index_cv);
    sum_ave2 = _global_av2(_index_cv);
    coutf << setw(10) << blk
          << setw(20) << setprecision(10) << average
          << setw(20) << setprecision(10) << sum_average/double(blk)
          << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if (_measure_chi){
    coutf.open("../OUTPUT/susceptibility.dat",ios::app);
    average  = _average(_index_chi);
    sum_average = _global_av(_index_chi);
    sum_ave2 = _global_av2(_index_chi);
    coutf << setw(10) << blk
          << setw(20) << setprecision(10) << average
          << setw(20) << setprecision(10) << sum_average/double(blk)
          << setw(20) << setprecision(10) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << setw(10) << blk << setw(20) << fraction << endl;
  coutf.close();
  
  return;
}

double System :: error(double acc, double acc2, int blk){ // Data blocking
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();
  return;
}

///////////////////////////////////////////////////////////////////

void System :: finalize(){ // Finalize the system
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}

void System :: write_configuration(){ // Write current configuration as a .xyz file in directory ../OUTPUT/CONFIG/
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open("../OUTPUT/CONFIG/config.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " 
              << setw(16) << setprecision(8) << _particle(i).getposition(0,true)/_side(0)          // x
              << setw(16) << setprecision(8) << _particle(i).getposition(1,true)/_side(1)          // y
              << setw(16) << setprecision(8) << _particle(i).getposition(2,true)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    this->write_velocities();
  } else {
    coutf.open("../OUTPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

void System :: write_velocities(){ // Write final velocities in an output file
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/velocities.out");
  if(coutf.is_open()){
    for(int i=0; i<_npart; i++){
      coutf << setw(16) << setprecision(8) << _particle(i).getvelocity(0)          // vx
            << setw(16) << setprecision(8) << _particle(i).getvelocity(1)          // vy
            << setw(16) << setprecision(8) << _particle(i).getvelocity(2) << endl; // vz
    }
  } else cerr << "PROBLEM: Unable to open velocities.dat" << endl;
  coutf.close();
  return;
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
