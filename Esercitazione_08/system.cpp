#include "system.h"

using namespace std;

///////////////////////////////////////////////////////////////////

int System :: get_nbl(){
    return _nblocks;
}

int System :: get_nsteps(){
    return _nsteps;
}

void System :: set_beta(double beta){
    _temperature = 1./beta;
    _beta = beta;
}

void System :: set_temperature(double T){
    _temperature = T;
    _beta = 1./T;
}

void System :: set_delta_ext(double delta_ext){
    _delta_ext = delta_ext; 
}  

///////////////////////////////////////////////////////////////////

void System :: initialize(bool sa, int blocks, int steps, double delta_int, double x0, vector<double> parameters){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

    this->average_reset();
    this->block_reset();

    _nblocks = blocks;
    _nsteps = steps;
    _delta_int = delta_int;
    _x0 = x0;
    _sa = sa;

    int p1, p2;  //Read from Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    int seed[4];  //Read the seed of the RNG
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed,p1,p2);

    if(sa==0){
        ofstream couta("acceptance_int.dat"); // Set the heading line in file acceptance_int.dat
        couta << "#    N_BLOCK:      ACCEPTANCE:" << endl;
        couta.close();

        ofstream couth("hamiltonian.dat");
        couth << "#      BLOCK:       ACTUAL_H:            H_AVE:             ERROR:" << endl;
        couth.close();
    }else if(sa==1){
        ofstream couta("acceptance_ext.dat"); // Set the heading line in file acceptance_ext.dat
        couta << "#    N_BLOCK:      ACCEPTANCE:" << endl;
        couta.close();

        ofstream coute("energy.dat");
        coute << "#      STEP:           E_AVE:             ERROR:" << endl;
        coute.close();
    }

    // first evaluation of the energy before starting simulated annealing
    this->data_blocking(_x0, parameters);

    _old_energy =  _global_av/double(this->get_nbl());

    return;
}

///////////////////////////////////////////////////////////////////

void System :: data_blocking(double x0, vector<double> parameters){ // Data blocking

    double x = x0;
    //ofstream couts("sqr_mod.dat");
    
    for(int i=0; i < this->get_nbl(); i++){ // Loop over blocks        
        for(int j=0; j < this->get_nsteps(); j++){ // Loop over steps in a block
            x = this->move(x, parameters); // Perform a MC step
            //couts << x[0] << endl; // Am I sampling right? Plot probability distribution
            this->hamiltonian(parameters, x);   // Compute integrating_function
        }
        this->averages(i+1);
        this->block_reset();
    }

    //couts.close();
}

double System :: move(double xold,  vector<double> parameters){ // Propose a MC move 
    
    _nattempts_int++;     // Update number of attempts performed on the system               
    
    double shift;         // Store the proposed translation
    double xnew;          // Translated position

    shift = _rnd.Rannyu(-1.0,1.0) * _delta_int;  // Uniform distribution in [-_delta;_delta)
    xnew = xold + shift;                         // Shift

    bool decision = false;
    double acceptance;

    double new_value = psi_2(parameters, xnew);
    double old_value = psi_2(parameters, xold);

    acceptance = new_value/old_value;

    if(acceptance >=1){
        decision = true;
    } else if(_rnd.Rannyu() < acceptance ){
        decision = true; //Metropolis acceptance step
    }
    
    if(decision){ //Metropolis acceptance evaluation
        _naccepted_int++;
        return xnew;
    } else return xold;  // If translation is rejected, restore the old configuration

}

void System :: hamiltonian(vector<double> parameters, double x){ // Compute H
    double sigma = parameters[0];
    double mu = parameters[1];
    double term1 = exp(-(x - mu) * (x - mu) / (2. * sigma * sigma)) * ((x - mu) * (x - mu) - sigma * sigma);
    double term2 = exp(-(mu + x) * (mu + x) / (2. * sigma * sigma)) * ((mu + x) * (mu + x) - sigma * sigma);
    double kinetic = -0.5*((term1 + term2)/ (pow(sigma, 4)))/sqrt(psi_2(parameters, x));
    double potential = pow(x, 4)-2.5*pow(x,2);

    _block_av += kinetic+potential;
    
    return;
}

double System :: psi_2(vector<double> parameters, double x){
    double sigma = parameters[0];
    double mu = parameters[1];
    return pow(exp(-pow((x-mu),2)/(2.*pow(sigma,2)))+exp(-pow((x+mu),2)/(2.*pow(sigma,2))),2);
}

///////////////////////////////////////////////////////////////////

void System :: averages(int blk){ // Block average

    ofstream coutf;

    _average     = _block_av / double(_nsteps);
    _global_av  += _average;
    _global_av2 += _average * _average; 

    if(_sa == 0){
        coutf.open("hamiltonian.dat",ios::app);
        coutf << setw(10) << blk
          << setw(20) << setprecision(10) << _average
          << setw(20) << setprecision(10) << _global_av/double(blk)
          << setw(20) << setprecision(10) << this->error(_global_av, _global_av2, blk) << endl;
        coutf.close();

        double fraction_int;
        coutf.open("acceptance_int.dat",ios::app);
        if(_nattempts_int > 0) fraction_int = double(_naccepted_int)/double(_nattempts_int);
        else fraction_int = 0.0; 
        coutf << setw(10) << blk << setw(20) << fraction_int << endl;
        coutf.close();
    }

    return;
}

double System :: error(double acc, double acc2, int blk){ // Compute error
    if(blk <= 1) return 0.0;
    else return sqrt(fabs(acc2/double(blk) - pow(acc/double(blk),2))/double(blk));
}

void System :: block_reset(){ // Reset block accumulators to zero
    _block_av = 0.0;
    return;
}

void System :: average_reset(){ // Reset block accumulators to zero
    _average    = 0.0;
    _global_av  = 0.0;
    _global_av2 = 0.0;
    _naccepted_int = 0;
    _nattempts_int = 0;
    _naccepted_ext = 0;
    _nattempts_ext = 0;
    return;
}

///////////////////////////////////////////////////////////////////

vector<double> System :: move_parameters(vector<double> old_parameters, int step){ // Propose a MC move 
    
    _nattempts_ext++;     // Update number of attempts performed on the system               
    
    vector<double> shift;               // Store the proposed translation
    vector<double> new_parameters;      // Translated position
    vector<double> result;

    //cout << "old energy: " << _old_energy << endl;

    double old_v = boltzmann(_old_energy);

    double delta = _delta_ext/_beta; // Change delta with temperature

    if(delta<0.01){
        delta = 0.01;
    }

    //cout << "delta: " << delta << endl;
    //cout << "old value: " << old_v << endl;

    for(int i=0; i<old_parameters.size(); i++){
        shift[i] = _rnd.Rannyu(-1.0,1.0) * delta;               // Uniform distribution in [-_delta;_delta)
        new_parameters.push_back(old_parameters[i] + shift[i]);       // Shift
    }

    this->data_blocking(_x0, new_parameters);

    //cout << "New parameters: " << new_parameters[0] << ' ' << new_parameters[1] << endl;

    bool decision = false;
    double acceptance;

    //cout << "New energy: " << _global_av/double(this->get_nbl()) << endl;

    double new_v = boltzmann(_global_av/double(this->get_nbl()));
    
    acceptance = new_v/old_v;

    //cout << "acceptance: " << acceptance << endl;

    if(acceptance >=1){
        decision = true;
    } else if(_rnd.Rannyu() < acceptance ){
        decision = true; //Metropolis acceptance step
    }

    //cout << "Decision: " << decision << endl;
    
    if(decision){ //Metropolis acceptance evaluation
        _naccepted_ext++;
        _old_energy = _global_av/double(this->get_nbl());
        result = new_parameters;
    } else {result=old_parameters;}  // If translation is rejected, restore the old configuration

    double fraction_ext;
    ofstream coutf;
    coutf.open("acceptance_ext.dat",ios::app);
    if(_nattempts_ext > 0) fraction_ext = double(_naccepted_ext)/double(_nattempts_ext);
    else fraction_ext = 0.0; 
    coutf << setw(10) << step << setw(20) << fraction_ext << endl;
    coutf.close();

    coutf.open("energy.dat",ios::app);
        coutf << setw(10) << step
          << setw(20) << setprecision(10) << _global_av/double(this->get_nbl())
          << setw(20) << setprecision(10) << this->error(_global_av, _global_av2, this->get_nbl()) << endl;
    coutf.close();

    return result;
}

double System :: boltzmann(double energy){
    return exp(-_beta*energy);
}