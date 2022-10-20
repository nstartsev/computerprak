#include <stdio.h>
#include <mpi.h>
# include <fstream>
# include <iostream>
# include <string>
# include <cmath>

using namespace std;

int main ( int argc, char *argv[] );
double boundary_condition ( double x, double time );
double initial_condition ( double x, double time, int condition_type, double x_max);
double rhs ( double x, double time );
void timestamp ( );
void update ( int id, int p );

int main(int argc, char **argv) {

    int ierr, num_procs, my_id;
    double wtime;
    ofstream x_file;
    int problem_type = 0;
    /* find out MY process ID, and how many processes were started. */

    ierr = MPI_Init(&argc, &argv);

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (my_id == 0){
        cout << "ADVECTION MPI" << endl;
        cout << "solve the advection equation" << endl;
        cout << "Number of processors= " << num_procs;
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

    double cfl = 0.8;
    double *u;
    ofstream u_file;
    double *u_new;
    int i;
    int Nmax = 10;
    int n = Nmax/p;
    MPI_Status status;
    int tag;
    double time;

    double time_max = 20.0;
    double time_min = 0.0;
    double *x;
    double hx;
    ofstream x_file;
    double x_max = 1000.0;
    double x_min = -x_max;
    double a = 1.0;
    int condition_type = 2;
    int problem_type = 2;

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
//  Set the values of H at the initial time.
//
    time = time_min;
    u = new double[n+2];
    u_new = new double[n+2];
    for ( i = 0; i <= n+1; i++ )
    {
        u[i] = initial_condition ( x[i], time , condition_type, x_max);
    }
    hx = ( x_max - x_min ) / ( double ) ( p * n - 1 );
    double dt = hx*cfl/a;
//
//  Check the CFL condition, have processor 0 print out its value,
//  and quit if it is too large.
//


    if ( id == 0 )
    {
        cout << "\n";
        cout << "UPDATE\n";
        cout << "  CFL stability criterion value = " << cfl << "\n";;
        cout << " N = " << Nmax << endl;
    }

    if ( 1 < cfl )
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
    for ( i = 1; i <= n; i++ )
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
//  Update the velocity.
//
        if (problem_type == 0){
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
        }
        if (problem_type == 1)
        {
            for ( i = 1; i <= n; ++i){
                u_new[i] = u[i] - dt*(u[i]*(u[i]-u[i-1])/(hx));
            }
        }
        if (problem_type == 2)
        {
            for ( i = 1; i <= n; ++i){
                u_new[i] = 1./2*(u[i+1]+u[i-1]) - dt*a/(2*hx)*(u[i+1]-u[i-1]);
            }
        }


//
//
// static boundary condition
/*
        if ( 0 == id )
        {
            u_new[1] = boundary_condition ( x[1], dt);
        }
        if ( id == p - 1 )
        {
            u_new[n] = boundary_condition ( x[n], dt);
        }
*/
        //
// periodic boundary condition

//  Send H[0] from ID 0 to ID p-1.

//
        if ( id == 0)
        {
            tag = 2;
            MPI_Send ( &u_new[2], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD );
        }
//
//  Receive H[0] from ID 0 .
//
        if ( id == p-1 )
        {
            tag = 2;
            MPI_Recv ( &u_new[n+1], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status );
        }
        //  Send H[N] to ID 0.
//
        if ( id == p - 1 )
        {
            tag = 2;
            MPI_Send ( &u_new[n-1], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );
        }
//
//  Receive H[0] from ID p-1.
//
        if ( id == 0 )
        {
            tag = 2;
            MPI_Recv ( &u_new[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status );
        }

//

//
//  Update time.
//
        time = time + dt;

        for ( i = 0; i <= n+1; i++ )
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
    if (p == 1)
    {
        double diff = 0.0;
        double du;
        ofstream data_file;
        ofstream data_table_csv;
        data_file.open("data.dat",ofstream ::out | ofstream::app);
        data_table_csv.open("data_lf.csv",ofstream ::out | ofstream::app);
        data_file << "N=  " << Nmax << endl;

        data_table_csv << "N;";
        data_table_csv << "DIFF_CONTINIOUS;";
        data_table_csv << "DIFF_INTEGRAL;";
        data_table_csv << "DIFF_LINF;";
        data_table_csv << "DIFF_L2;";
        data_table_csv << "DIFF_L2/N" << endl;
        data_table_csv << Nmax << ';';
        //continious
        for ( i = 1; i <=n; i++)
        {
            diff += abs(u[i] - initial_condition(x[i], time, condition_type, x_max)) ;
        }
        data_file << "DIFF_CONTINIOUS=  " << diff << endl;
        data_table_csv << diff << ';';
        //integral
        diff = 0.0;
        for ( i = 1; i <=n; i++)
        {
            diff += abs(u[i] - initial_condition(x[i], time, condition_type, x_max))*hx;
        }
        data_file << "DIFF_INTEGRAL=  " << diff << endl;
        data_table_csv << diff << ';';
        //l_inf
        diff = 0.0;
        for ( i = 1; i <=n; i++)
        {
            diff = max(abs(u[i] - initial_condition(x[i], time, condition_type, x_max)),diff);
        }
        data_file << "DIFF_LINF=  " << diff << endl;
        data_table_csv << diff << ';';
        // l_2
        diff = 0.0;
        for ( i = 1; i <=n; i++)
        {
            du = abs(u[i] - initial_condition(x[i], time, condition_type, x_max));
            diff += du*du;
        }
        diff = pow(diff, 0.5);
        data_file << "DIFF_L2=  " << diff << endl;
        data_table_csv << diff << ';';
        // l_2/N
        diff = 0.0;
        for ( i = 1; i <=n; i++)
        {
            du = abs(u[i] - initial_condition(x[i], time, condition_type, x_max));
            diff += du*du;
        }
        diff = pow(diff/Nmax, 0.5);
        data_file << "DIFF_L2/N=   " << diff << endl;
        data_table_csv << diff << endl;
//        data_file.close();
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
double initial_condition ( double x, double time, int condition_type , double x_max)
{
    double value = 0.0;
    if (condition_type == 0)
    {
        if (abs(x) <= x_max/2){
            value = (x_max/2-abs(x));
        }
        else {
            value = 0.0;
        }
    }
    if (condition_type == 1){
        value = (1-abs(x)/x_max);
    }
    if (condition_type == 2)
    {
        value = cos(2*M_PI/(2*x_max)*x);
    }
    return value;
}
double boundary_condition ( double x, double time )

{
    double value;
    value = 0.0;
    return value;
}
double rhs ( double x, double time )

{
    double value;
    value = 0.0;
    return value;
}