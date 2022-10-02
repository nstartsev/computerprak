#include <stdio.h>
#include <mpi.h>
# include <fstream>
# include <iostream>
# include <string>
# include <cmath>

using namespace std;

int main ( int argc, char *argv[] );
double boundary_condition ( double x, double time );
double initial_condition ( double x, double time );
double rhs ( double x, double time );
void timestamp ( );
void update ( int id, int p );

int main(int argc, char **argv) {

    int ierr, num_procs, my_id;
    double wtime;
    ofstream x_file;
    /* find out MY process ID, and how many processes were started. */

    ierr = MPI_Init(&argc, &argv);

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (my_id == 0){
        cout << "ADVECTION MPI" << endl;
        cout << "solve the advection equation" << endl;
    }
    if ( my_id == 0 )
    {
        wtime = MPI_Wtime ( );
    }

    update(my_id, num_procs);

    if ( my_id == 0 )
    {
        wtime = MPI_Wtime ( ) - wtime;

        cout << "\n";
        cout << "  Wall clock elapsed seconds = " << wtime << "\n";
    }
    ierr = MPI_Finalize();

}
void update(int id, int p){

    double cfl;
    double *u;
    ofstream u_file;
    double *u_new;
    int i;
    int n = 33;
    MPI_Status status;
    int tag;
    double time;
    double dt = 0.05;
    double time_max = 10.0;
    double time_min = 0.0;
    double *x;
    double hx;
    ofstream x_file;
    double x_max = 16.0;
    double x_min = -16.0;
    double a = 1.0;

//
//  Set the X coordinates of the N nodes.
//  We don't actually need ghost values of X but we'll throw them in
//  as X[0] and X[N+1].
//
    x = new double[n+2];

    for ( i = 0; i <= n + 1; i++ )
    {
        x[i] = ( ( double ) (         id * n + i - 1 ) * x_max
                 + ( double ) ( p * n - id * n - i     ) * x_min )
               / ( double ) ( p * n              - 1 );
    }
//
//
//  Set the values of H at the initial time.
//
    time = time_min;
    u = new double[n+2];
    u_new = new double[n+2];
    for ( i = 0; i <= n+1; i++ )
    {
        u[i] = initial_condition ( x[i], time );
    }
    hx = ( x_max - x_min ) / ( double ) ( p * n - 1 );
//
//  Check the CFL condition, have processor 0 print out its value,
//  and quit if it is too large.
//
    cfl = a*dt/hx;

    if ( id == 0 )
    {
        cout << "\n";
        cout << "UPDATE\n";
        cout << "  CFL stability criterion value = " << cfl << "\n";;
    }

    if ( 1 <= cfl )
    {
        if ( id == 0 )
        {
            cout << "\n";
            cout << "UPDATE - Warning!\n";
            cout << "  Computation cancelled!\n";
            cout << "  CFL condition failed.\n";
            cout << "  1 <= a*dt/hx =  " << cfl << "\n";
        }
        return;
    }
//
//  Write out the values of H.
//
    u_file.open ( "u_data" + to_string(id) + ".dat" );
    u_file << "TITLE= 'SOLUTION'" << endl;
    u_file << "VARIABLES= " << "'X'" << '\t' << "'U'" << endl;
    u_file << "ZONE T= " << "'" << id << "'" << '\t' << "I= " << n << '\t' << "F=POINT" << endl;
    u_file << "SOLUTIONTIME= " << time << endl;
    for ( i = 0; i <= n+1; i++ )
    {
        u_file << setprecision (3) << fixed << x[i] <<  '\t' << u[i] << endl;
    }
    u_file << "\n";
//
//  Compute the values of H at the next time, based on current data.
//
    while (time < time_max)
    {

        if (id == 1)
        {
//            cout << time << '\t' << id << '\t' << u[0]  << endl;
        }

//
//  Send H[1] to ID-1.
//
        if ( 0 < id )
        {
            tag = 1;
            MPI_Send ( &u[1], 1, MPI_DOUBLE, id-1, tag, MPI_COMM_WORLD );
        }
//
//  Receive H[N+1] from ID+1.
//
        if ( id < p-1 )
        {
            tag = 1;
            MPI_Recv ( &u[n+1], 1,  MPI_DOUBLE, id+1, tag, MPI_COMM_WORLD, &status );
        }

//
//  Send H[N] to ID+1.
//
        if ( id < p-1 )
        {
            tag = 2;
            MPI_Send ( &u[n], 1, MPI_DOUBLE, id+1, tag, MPI_COMM_WORLD );
        }
//
//  Receive H[0] from ID-1.
//
        if ( 0 < id )
        {
            tag = 2;
            MPI_Recv ( &u[0], 1, MPI_DOUBLE, id-1, tag, MPI_COMM_WORLD, &status );
        }


//
//  Update the temperature based on the four point stencil.
//
        if (a>0){
            for ( i = 1; i <= n; ++i){
                u_new[i] = u[i] - dt*(a*(u[i]-u[i-1])/(hx));
            }
        }
        else{
            for ( i = 1; i <= n; ++i){
                u_new[i] = u[i] - dt*(a*(u[i+1]-u[i])/(hx));
            }
        }

//
//  H at the extreme left and right boundaries was incorrectly computed
//  using the differential equation.  Replace that calculation by
//  the boundary conditions.
//
        if ( 0 == id )
        {
            u_new[1] = boundary_condition ( x[1], dt);
        }
        if ( id == p - 1 )
        {
            u_new[n] = boundary_condition ( x[n], dt);
        }
//
//  Update time and temperature.
//

        time = time + dt;

        for ( i = 1; i <= n; i++ )
        {
            u[i] = u_new[i];
        }

//
//  In single processor mode, add current solution data to output file.
//

        u_file << "ZONE T= " << "'" << id << "'" << '\t' << "I= " << n << '\t' << "F=POINT" << endl;
        u_file << "SOLUTIONTIME= " << time << endl;
        for ( i = 1; i <= n; i++ )
        {
            u_file << setprecision (3) << fixed << x[i] <<  '\t' << u[i] << endl;
        }
        u_file << "\n";

    }

    if ( p == 1 )
    {
        u_file.close ( );
    }

    delete [] u;
    delete [] u_new;
    delete [] x;

    return;
}
double initial_condition ( double x, double time )

{
    double value;

    if (abs(x) <= 1){
        value = (1-abs(x));
    }
    else {
        value = 0.0;
    }

    return value;
}
double boundary_condition ( double x, double time )

{
    double value;
//
//  Left condition:
//
    value = 0.0;
    return value;
}
double rhs ( double x, double time )

{
    double value;

    value = 0.0;

    return value;
}