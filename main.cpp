#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

int main() {
    //mesh
    double Xmax = 10;
    int Im = 500;
    int Imax = 1000;
    vector <double> x(Imax,0);
    double hx = 2*Xmax/(Im);
    for (size_t i = 0; i <= Im; ++i){
        x[i] = hx*i - Xmax;
    }

    // init u0(x)
    vector <double> u0(Imax,0);
    for (size_t i = 0; i <= Im; ++i){
        if (abs(x[i]) <= 1){
            u0[i] = 1-x[i]*x[i];
        }
    }

    vector <double> u(Imax, 0);
    for (size_t i = 0 ; i <= Im; ++i){
        u[i] = u0[i];
    }

    double Tmax = 10;
    double dt = 0.05;

}
