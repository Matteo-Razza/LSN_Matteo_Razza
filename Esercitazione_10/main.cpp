#include <iostream>
#include "city.h"
#include "path.h"
#include "salesman.h"
#include <mpi.h>

using namespace std;

int main (int argc, char *argv[]){

    int n_cities = 110;
    vector<city> cities(n_cities);
    int population = 3000;
    int generations = 4000;
    int migration = 40;

    const int itag = 1;
    MPI_Status stat;

    // Reading positions from a file
    double x;
    double y;
    
    ifstream filein("cap_prov_ita.dat");
    for(int j=0; j<n_cities; j++){
        filein >> x >> y;
        cities[j].initialize(x, y);
    }
    filein.close();

    // Start MPI
    int size, rank;
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //cout<<" Sono il nodo "<<rank<<" dei "<< size <<" che hai utilizzato!"<<endl;

    Random rnd; 
    int p1, p2;  // Read from Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("Primes");
    for (int i=0; i< rank+1; ++i){
        Primes >> p1 >> p2; //This is crucial, different values for each rank
    }
    Primes.close();
    int seed[4];  // Read the seed of the RNG
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);

    salesman TS;

    // Creation of the starting population
    TS.initialize(cities, rnd, population);

    ofstream fastest_out("fastest_rank_"+ std::to_string(rank) +".dat");
    ofstream path_out("path_rank_"+std::to_string(rank)+".dat");

    for(int i=0; i<generations-1; i++){
        // Sort path
        TS.sort_paths();

        fastest_out << i+1 << " " << TS.cout_length_single_path(0) << endl;

        // Writing positions on a file every once in a while
        if((i+1)%(generations/10)==0){
            for(int j=0; j<n_cities; j++){
                path best_path = TS.get_best_path();
                vector<int> best = best_path.get_path();
                path_out << j+1 << " " << setprecision(18) << cities[best[j]-1].get_y() << " " << setprecision(18) << cities[best[j]-1].get_x() << endl;
                if(j==n_cities-1){
                    path_out << n_cities+1 << " " << setprecision(18) << cities[best[0]-1].get_y() << " " << setprecision(18) << cities[best[0]-1].get_x() << endl;
                }
            }
        }

        if((i+1)%migration==0){

            int europe = 0;
            int america = 0;
            path best_eu;
            path best_am;
            vector<int> best_europe(n_cities);
            vector<int> best_america(n_cities);
            path recived;

            if (rank == 0) {
                cout << "Migration number " << (i+1)/migration << endl;
                europe = int(rnd.Rannyu(0, size)); 
                america = int(rnd.Rannyu(0, size)); 
                while (america == europe){
                    america = int(rnd.Rannyu(0, size));
                }
            }
            
            // Broadcast: everyone knows who is switching
            MPI_Bcast(&europe, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&america, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            if(rank == europe){
                best_eu = TS.get_best_path();
                best_europe = best_eu.get_path();
                MPI_Send(best_europe.data(), best_europe.size(), MPI_INT, america, itag, MPI_COMM_WORLD); 
                MPI_Recv(best_america.data(), best_america.size(), MPI_INT, america, itag, MPI_COMM_WORLD, &stat);
                recived.initialize(best_america, n_cities);
                TS.set_path(recived, population-1);
            }

            if(rank == america){
                best_am = TS.get_best_path();
                best_america = best_am.get_path();
                MPI_Recv(best_europe.data(), best_europe.size(), MPI_INT, europe, itag, MPI_COMM_WORLD, &stat);
                MPI_Send(best_america.data(), best_america.size(), MPI_INT, europe, itag, MPI_COMM_WORLD);
                recived.initialize(best_europe, n_cities);
                TS.set_path(recived, population-1);
            }

            if(rank == europe or rank == america){
                cout << "Rank: " << rank << ", lenght first path (0): " << TS.cout_length_single_path(0) << endl;
                cout << "Rank: " << rank << ", lenght last path (population-1): " << TS.cout_length_single_path(population-1) << endl;
            }

            // Send-Recive-->Recive-Send: should be safe
        }

        vector<path> result(population);

        for(int j=0; j<population; j=j+2){
            // Select parents
            int selected_father = TS.select_parents();
            int selected_mother = TS.select_parents();

            if(rnd.Rannyu()<0.8){ // Crossover al 80%
                // Generate children
                path son = TS.crossover(selected_father, selected_mother);
                path daughter = TS.crossover(selected_mother, selected_father);

                result[j]=son;
                result[j+1]=daughter;
            } else{
                // Keep parents
                result[j]=TS.get_path(selected_father);
                result[j+1]=TS.get_path(selected_mother);
            }
        }

        TS.set_paths(result);

        for(int j=0; j<population; j++){
            // Mutations
            if(rnd.Rannyu()<0.1){
                TS.shift(j);
            }
            if(rnd.Rannyu()<0.1){
                TS.permutate(j);
            }
            if(rnd.Rannyu()<0.2){
                TS.invert(j);
            }
            if(rnd.Rannyu()<0.2){
                TS.swap(j);
            }
        }
    }

    TS.sort_paths();

    cout << "Rank: " << rank << ", lenght best path: " << TS.cout_length_single_path(0) << endl;
    cout << "Rank: " << rank << ". ";
    TS.cout_single_path(0);

    fastest_out << population << " " << TS.cout_length_single_path(0) << endl;

    for(int j=0; j<n_cities; j++){
        path best_path = TS.get_best_path();
        vector<int> best = best_path.get_path();
        path_out << j+1 << " " << setprecision(18) << cities[best[j]-1].get_y() << " " << setprecision(18) << cities[best[j]-1].get_x() << endl;
        if(j==n_cities-1){
            path_out << n_cities+1 << " " << setprecision(18) << cities[best[0]-1].get_y() << " " << setprecision(18) << cities[best[0]-1].get_x() << endl;
        }
    }

    path_out.close();
    fastest_out.close();

    MPI_Finalize();

    return 0;
}

// mpirun -np 2 ./mpi_main.exe