#ifndef __salesman__
#define __salesman__

#include "path.h"
#include <algorithm>

using namespace std;

class salesman {

private:
    vector<city> _map;
    vector<path> _itinera;
    Random _rnd;

public:
    
    void initialize(const vector<city> cities, Random random, int population){
        _map = cities;
        _rnd = random;
        int l=0;
        int k=0;

        _itinera.resize(population);
        vector<int> vec(_map.size());

        for(int i=0; i<population; i++){
            
            for(int a=0; a<_map.size(); a++){
                vec[a]=a+1;
            }

            for(int j=0; j<17; j++){
                l = int(_rnd.Rannyu(1,_map.size())); // extract l between 1 and _map.size()-1
                k = int(_rnd.Rannyu(1,_map.size())); // extract k between 1 and _map.size()-1
                std::swap(vec[l], vec[k]);
            }

            _itinera[i].initialize(vec, _map.size());

            if(_itinera[i].check()==false){
                cerr << "Starting population don't fulfil the bonds. Exit." << endl;
                exit(-1);
            }
        }
    }

    int select(){
        int selected = int(_itinera.size()*pow(_rnd.Rannyu(),0.5));

        if(selected<0 or selected >_itinera.size()-1){
            cout << "Selection failed" << endl;
            exit(-1);
        }

        return selected;
    }

    int select_parents(){
        int selected = int(_itinera.size()*pow(_rnd.Rannyu(), 1.5));

        if(selected<0 or selected >_itinera.size()-1){
            cout << "Selection failed" << endl;
            exit(-1);
        }

        return selected;
    }

    void swap(int selected){
        int l = int(_rnd.Rannyu(1,_map.size())); // extract l between 1 and _map.size()-1
        int k = int(_rnd.Rannyu(1,_map.size())); // extract k between 1 and _map.size()-1
        //cout << "l = " << l << ", k = " << k << endl; 
        std::swap(_itinera[selected].get_path()[l], _itinera[selected].get_path()[k]);
    }

    void shift(int selected){
        vector<int> vec(_map.size());
        vec = _itinera[selected].get_path();

        int i = int(_rnd.Rannyu(1,_map.size())); //maybe greater
        int m = int(_rnd.Rannyu(1,_map.size()-1));
        int n = int(_rnd.Rannyu(1,_map.size()-m));

        //cout << "i = " << i << ", m = " << m << ", n = " << n << endl; 

        vector<int> copy(_map.size()-1);
        for(int k=0; k<_map.size(); k++){
            if(i+k<_map.size()){
                copy[k]=vec[i+k];
            }else if(i+k>_map.size()){
                copy[k-1]=vec[i+k-_map.size()];
            }
        }

        rotate(copy.begin(), copy.begin()+m, copy.begin()+m+n);

        vec[0]=1;
        for(int k=0; k<_map.size(); k++){
            if(i+k<_map.size()){
                vec[i+k]=copy[k];
            }else if(i+k>_map.size()){
                vec[i+k-_map.size()]=copy[k-1];
            }
        }       

        _itinera[selected].set_path(vec);

        if(_itinera[selected].check()==false){
            cout << "Problem during shift" << endl;
            exit(-1);
        }
    }

    void permutate(int selected){
        vector<int> vec(_map.size());
        vec = _itinera[selected].get_path();

        int m = int(_rnd.Rannyu(1,17));

        int i = int(_rnd.Rannyu(1, _map.size()));
        int j = int(_rnd.Rannyu(i+m, i+_map.size()-1-m));

        vector<int> copy(_map.size()-1);
        for(int k=0; k<_map.size(); k++){
            if(i+k<_map.size()){
                copy[k]=vec[i+k];
            }else if(i+k>_map.size()){
                copy[k-1]=vec[i+k-_map.size()];
            }
        }

        swap_ranges(copy.begin(), copy.begin()+m, copy.begin()+(j-i));

        vec[0]=1;
        for(int k=0; k<_map.size(); k++){
            if(i+k<_map.size()){
                vec[i+k]=copy[k];
            }else if(i+k>_map.size()){
                vec[i+k-_map.size()]=copy[k-1];
            }
        }

        //cout << "i = " << i << ", m = " << m << ", j = " << j-(_map.size()-1) << endl;        

        _itinera[selected].set_path(vec);

        if(_itinera[selected].check()==false){
            cout << "Problem during permutation" << endl;
            exit(-1);
        }
    }

    void invert(int selected){ 
        vector<int> vec(_map.size());
        vec = _itinera[selected].get_path();
        int m = int(_rnd.Rannyu(2, _map.size()));
        int i = int(_rnd.Rannyu(1, _map.size()));
        vector<int> copy(m);
        for(int j=i; j<i+m; j++){
            if(j<_map.size()){
                copy[j-i]=vec[j];
            }
            else if(j>=_map.size()){
                copy[j-i]=vec[j-(_map.size()-1)];
            }
        }
        reverse(copy.begin(), copy.end());
        for(int j=i; j<i+m; j++){
            if(j<_map.size()){
                vec[j]=copy[j-i];
            }
            else if(j>=_map.size()){
                vec[j-(_map.size()-1)]=copy[j-i];
            }
        }
        //cout << "i = " << i << ", m = " << m << endl;        

        _itinera[selected].set_path(vec);

        if(_itinera[selected].check()==false){
            cout << "Problem during invertion" << endl;
            exit(-1);
        }
    }

    void sort_paths() {
        sort(_itinera.begin(), _itinera.end(), [&](const path& a, const path& b) {
            return a.length(_map) < b.length(_map);
        });
    }

    path crossover(int father, int mother){

        //cout << father << " " << mother << endl;

        int divider = _rnd.Rannyu(1,_map.size());
        vector<int> son(_map.size());
        //vector<int> daughter(_map.size());

        vector<int> dad(_map.size());
        dad = _itinera[father].get_path();

        vector<int> mom(_map.size());
        mom = _itinera[mother].get_path();

        int index_son=divider;
        //int index_daughter=divider;

        son[0] = 1;
        //daughter[0] = 1;

        unordered_set<int> included_son;
        included_son.insert(1);
        //unordered_set<int> included_daughter;
        //included_daughter.insert(1);

        // Copy the first part from one parent
        for (int i = 1; i < divider; i++) {
            son[i] = dad[i];
            included_son.insert(dad[i]);
            //daughter[i] = mom[i];
            //included_daughter.insert(mom[i]);
        }
    
        // Fill in the rest of the cities from the consort parent in order
        for (int i=0; i<_map.size(); i++) {
            if (included_son.find(mom[i])==included_son.end()){ // If element not found, included_son.find(...) returns included.end()
                son[index_son] = mom[i];
                index_son ++;
            }

            // if (included_daughter.find(dad[i])==included.end()){ // If element not found, included_daughter.find(...) returns included.end()
            //     son[index_daughter++] = dad[i];
            // }
        }

        path result;

        result.initialize(son, _map.size());

        if(result.check() == true){
            return result;
        }else{
            cout << "Problem during crossover operations" << endl;
            exit(-1);
        }
    }

    void new_generation(path son, path daughter, int i, int j){
        _itinera[i]=son;
        _itinera[j]=daughter;
        if(_itinera[i].check()==false or _itinera[j].check()==false){
            cout << "Problem during creation of new generation" << endl;
            exit(-1);
        }
    }

    void cout_length_paths() const{
        for(int i=0; i<_itinera.size(); i++){
            cout << "L_2(x" << i+1 << "): " << _itinera[i].length(_map) << endl;
        }
    }

    double cout_length_single_path(int i) const{
        return _itinera[i].length(_map);
    }

    void cout_single_path(int i) const{
        _itinera[i].cout_path();
    }

    const vector<path>& get_paths() const{
        return _itinera;
    }

    vector<path>& get_paths(){
        return _itinera;
    }
    
    path get_best_path(){
        return _itinera[0];
    }

    path get_path(int i){
        return _itinera[i];
    }

    void set_paths(vector<path> vec){
        _itinera = vec;
    }

};

#endif // __salesman__