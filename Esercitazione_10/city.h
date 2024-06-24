#ifndef __city__
#define __city__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;

class city {
private:
    double _x;
    double _y;
public: 
    void initialize(double x, double y){
        _x = x;
        _y = y;
    }
    double distance(const city& one, const city& two) const{
        double distance = sqrt(pow(one.get_x()-two.get_x(), 2)+pow(one.get_y()-two.get_y(), 2));
        return distance;
    }
    double get_x() const{
        return _x;
    }
    double get_y() const{
        return _y;
    }
    void set_x(double x){
        _x = x;
    }
    void set_y(double y){
        _y = y;
    }

};

#endif // __city__
