#ifndef __path__
#define __path__

#include "city.h"
#include "random.h"
#include <unordered_set>

using namespace std;

class path {

private:
    vector<int> _path;
    int _dim; 

public: 
    
    void initialize(const vector<int> vec, int dim){
        
        for(int i=0; i<dim; i++){
            _path.push_back(vec[i]);
        }
    
        _dim = dim;
    }
    
    bool check() {
        if (_path.size() != _dim){
            cout << "Dimension problem" << endl;
            return false;
        }

        std::unordered_set<int> labels;
        int stop;

        if(_path[0]!=1){
            cout << "First element not 1" << endl;
            return false;
        }

        for(int i=1; i<_dim; i++){
            stop = _path[i];
            if(!labels.insert(stop).second) {
                cout << "Twice same stop" << endl;
                return false;
            }
        }
    
        return labels.size() == _dim-1;
    }
   
    double length(const vector<city>& cities) const{
        
        double length = 0.0;

        for(int i=0; i<_dim-1; i++){
            length += sqrt(pow(cities[_path[i]-1].get_x()-cities[_path[i+1]-1].get_x(), 2)+pow(cities[_path[i]-1].get_y()-cities[_path[i+1]-1].get_y(), 2));
        }

        length += sqrt(pow(cities[_path[_dim-1]-1].get_x()-cities[_path[0]-1].get_x(), 2)+pow(cities[_path[_dim-1]-1].get_y()-cities[_path[0]-1].get_y(), 2));

        return length;
    }

    vector<int>& get_path(){
        return _path;
    }

    void set_path(vector<int> iter){
        _path = iter;
    }

    void cout_path() const{
        cout << "Path: ";
        for(int i=0; i<_path.size(); i++){
            cout << _path[i] << " ";
        }
        cout << endl;
    }

};

#endif // __path__
