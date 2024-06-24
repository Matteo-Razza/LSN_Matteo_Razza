#include <iostream>
#include "city.h"
#include "path.h"
#include "salesman.h"

using namespace std;

int main (int argc, char *argv[]){

    Random rnd; 
    int p1, p2;  // Read from Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    int seed[4];  // Read the seed of the RNG
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);

    salesman TS;
    vector<city> cities(34);
    int population = 3000;
    int generations = 500;

    // // Placing 34 cities on a circumference of radius 1
    // double theta;
    // for(int i=0; i<34; i++){
    //     theta = rnd.Rannyu(0, 2.*M_PI);
    //     cities[i].initialize(cos(theta), sin(theta));
    // }

    // Placing 34 cities in a square
    double x;
    double y;
    for(int i=0; i<34; i++){
        x = rnd.Rannyu(-1, 1);
        y = rnd.Rannyu(-1, 1);
        cities[i].initialize(x, y);
    }

    // Writing positions on a file
    ofstream out("squ_map.dat");
    for(int j=0; j<34; j++){
        out << j+1 << " " << cities[j].get_x() << " " << cities[j].get_y() << endl;
    }
    out.close();

    // Creation of the starting population
    TS.initialize(cities, rnd, population);

    //ofstream fastest_out("squ_fastest.dat");
    ofstream average_out("squ_average.dat");

    for(int i=0; i<generations; i++){

        // Sort path
        TS.sort_paths();
        cout << TS.cout_length_single_path(0) << endl;

        //fastest_out << i+1 << " " << TS.cout_length_single_path(0) << endl;

        double average = 0;
        for(int p=0; p<population/2; p++){
            average += TS.cout_length_single_path(p);
        }

        average_out << i+1 << " " << average/(double(population)/2.) << endl;

        // // Writing positions on a file
        // ofstream path_out("Results/Square/path_"+std::to_string(i+1));
        // for(int j=0; j<34; j++){
        //     path best_path = TS.get_best_path();
        //     vector<int> best = best_path.get_path();
        //     path_out << j+1 << " " << cities[best[j]-1].get_x() << " " << cities[best[j]-1].get_y() << endl;
        //     if(j==33){
        //         path_out << 35 << " " << cities[best[0]-1].get_x() << " " << cities[best[0]-1].get_y() << endl;
        //     }
        // }
        // path_out.close();

        vector<path> result(population);

        for(int j=0; j<population; j=j+2){
            // Select parents
            int selected_father = TS.select_parents();
            int selected_mother = TS.select_parents();

            if(rnd.Rannyu()<0.7){ // Crossover al 70%
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
            if(rnd.Rannyu()<0.15){
                TS.shift(j);
            }
            if(rnd.Rannyu()<0.15){
                TS.permutate(j);
            }
            if(rnd.Rannyu()<0.15){
                TS.invert(j);
            }
            if(rnd.Rannyu()<0.15){
                TS.swap(j);
            }
        }
    }

    TS.sort_paths();
    cout << TS.cout_length_single_path(0) << endl;
    TS.cout_single_path(0);

    // for(int k=0; k<200; k++){
    //     TS.cout_single_path(k);
    // }

    //fastest_out.close();
    average_out.close();

    return 0;
}