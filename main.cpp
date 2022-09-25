#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

int main() {
    //mesh
    double Xmax = 20;
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
            u0[i] = (1-abs(x[i]));
        }
    }

    vector <double> u(Imax, 0);
    vector <double> u_new(Imax, 0);
    for (size_t i = 0 ; i <= Im; ++i){
        u[i] = u0[i];
    }
    // time and problem parameters
    double Tmax = 10;
    double dt = 0.05;
    double t = 0;
    double a = 1.;

    ofstream sol;
    sol.open("solution.dat");
    sol << "TITLE= 'SOLUTION'" << endl;
    sol << "VARIABLES= " << "'X'" << '\t' << "'U'" << endl;

    cout << hx << endl << a*dt << endl;
    // algorithm
    while (t <= 10.0){
        //print solution
        sol << "ZONE T= " << "'" << t << "'" << '\t' << "I= " << Im << '\t' << "F=POINT" << endl;
        sol << "SOLUTIONTIME= " << t << endl;
        for (size_t i = 0 ; i < Im; ++i){
            sol << x[i] << '\t' << u[i] << endl;
        }
        for (size_t i = 1; i < Im; ++i){
            u_new[i] = u[i] - dt*(a*(u[i]-u[i-1])/(hx));
        }
        for (size_t i = 1; i < Im; ++i){
            u[i] = u_new[i];
        }
        t += dt;

    }
    sol.close();
}
