#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main(int argc, char **argv) {

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
    int initial_condition = 1;
    vector <double> u0(Imax,0);

    if (initial_condition == 0){
        for (size_t i = 0; i <= Im; ++i){
            if (abs(x[i]) <= 1){
                u0[i] = (1-abs(x[i]));
            }
        }
    }
    if (initial_condition == 1){
        for (size_t i = 0; i <= Im; ++i){
                u0[i] = (1-abs(x[i])/Xmax);
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
    int problem = 1;

    ofstream sol;
    sol.open("solution.dat");
    sol << "TITLE= 'SOLUTION'" << endl;
    sol << "VARIABLES= " << "'X'" << '\t' << "'U'" << endl;

    cout << "hx= " << hx << endl << "adt= " <<  abs(a*dt) << endl;
    // algorithm
    while (t <= Tmax){
        //print solution
        sol << "ZONE T= " << "'" << t << "'" << '\t' << "I= " << Im << '\t' << "F=POINT" << endl;
        sol << "SOLUTIONTIME= " << t << endl;
        for (size_t i = 0 ; i < Im; ++i){
            sol << x[i] << '\t' << u[i] << endl;
        }
        if (problem == 0){
            if (a>0){
                for (size_t i = 1; i < Im; ++i){
                    u_new[i] = u[i] - dt*(a*(u[i]-u[i-1])/(hx));
                }
            }
            else{
                for (size_t i = 1; i < Im; ++i){
                    u_new[i] = u[i] - dt*(a*(u[i+1]-u[i])/(hx));
                }
            }
        }
        if (problem == 1){
            if (a>0){
                for (size_t i = 1; i < Im; ++i){
                    u_new[i] = u[i] - dt*(u[i]*(u[i]-u[i-1])/(hx));
                }
            }
            else{
                for (size_t i = 1; i < Im; ++i){
                    u_new[i] = u[i] - dt*(u[i]*(u[i+1]-u[i])/(hx));
                }
            }
        }
        for (size_t i = 1; i < Im; ++i) {
            u[i] = u_new[i];
        }
        t += dt;

    }
    sol.close();
}
