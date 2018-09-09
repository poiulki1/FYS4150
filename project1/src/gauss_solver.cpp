#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "time.h"

#include <armadillo>

using namespace std;
using namespace arma;


inline double f(double x){return 100.0*exp(-10.0*x);} //function we want to solve  --RHS
inline double exact(double x){return 1.0 - (1 - exp(-10))*x - exp(-10*x);} //exact function --found analytically --needed for error calculation

double rel_error(double, double);
void general_solver(int, double*, double*, double*, double*, double*);
void specific_solver(int, double*, double*, double*);

int main(int argc, char* argv[]){
    int N, algo, exponent;
    if(argc < 3){
        cout << "You need to give two arguments: \n first: n - 10^n to calculate number of mesh points \n second: algo: general - 0, spec - 1" << endl;
        exit(1);
    }
    else{
        exponent = atoi(argv[1]);
        algo = atoi(argv[2]);
    }
    ofstream myfile_sim ("sim_data.txt"); //simulated data - depends on the second cml arg. which algo will be used to write the data
    ofstream myfile_error ("error_data.txt"); //error data - to plot relative error with python

    for(int i = 1; i <= exponent; i++){
        N = (int) pow(10.0,i);

    myfile_sim << N << endl; //writing this to sim data file so we can use this to slicing the data for each N when plotting with python
    double dt = 1.0/(N+1); //time step

    //allocating memory dynamicly for our vectors which will be used to calculate solution
    double* a = new double[N+1]; //sub main diagonal
    double* b = new double[N+2]; //main diagonal
    double* c = new double[N+1]; //above main diagonal
    double* s = new double[N+2]; //our RHS function --f(x)*h^2
    double* v = new double[N+2]; //simulated solution

    double* err = new double[N+1]; //error array

    //initializing the arrays/vectors
    for(int i = 0; i < N+2; i++){
        s[i] = dt*dt*f(dt*i);
        a[i] = -1.0;
        b[i] = 2.0;
        c[i] = -1.0;
    }
    //Dirichlet boundary conditions
    s[0]= s[N+2] = 0;

    //the last main diagonal element
    b[N+2] = 2.0;

    //declearing some time variables to calculate time used by different algorithms
    clock_t start_time, finish_time;
    double time_used;

    //algorithm choice - depends on cml arg.
    if(algo==0){
        start_time = clock();
        general_solver(N, a, b, c, s, v);
        finish_time = clock();
        time_used = ((double) (finish_time - start_time))/CLOCKS_PER_SEC;
        cout << "General algorithm uses: " << time_used << " seconds, with N= " << N << endl;
    }
    else if(algo==1){
        start_time = clock();
        specific_solver(N,b,s,v);
        finish_time = clock();
        time_used = ((double) (finish_time - start_time))/CLOCKS_PER_SEC;
        cout << "Specific algorithm uses: " << time_used << " seconds, with N= " << N << endl;
    }

    //else if(algo==2){
    //    LU_decomp
    //}

    //calculating the relative error and finding the biggest value - needed for later analysis/plot
    for(int k = 1; k < N+2; k++){

        err[k] = rel_error(v[k], exact(dt*k));
     }
        double biggest_error = err[0]; //starting with the first error as the biggest one

     for(int j = 1; j < N; j++){
        if(abs(err[j]) >= abs(biggest_error)){
                biggest_error = err[j];
         }
      }
     //writing to file
      myfile_error << setprecision(12) << left << setw(25) << biggest_error
                                << left << setw(25) << N
                                << endl;
      //writing all the values of the simulated solution with corresponding mesh-points to the file (for all of the N)
      for(int p = 0; p <= N+1; p++){
      myfile_sim << setprecision(8) << left << setw(20) << v[p]
                                << left << setw(20) << p*dt
                                << endl;}


    }
}

//general gaussian elimination - all the elements can have different different values
void general_solver(int N, double* a, double* b, double* c, double* s, double* v){
    //forward
    for(int i = 2; i <= N+1; i++){
        b[i] -= (a[i-1]*c[i-1])/b[i-1];
        s[i] -= s[i-1]*(a[i-1]/b[i-1]);

    }
    //backward
    v[N-1] = s[N-1]/b[N-1];

    for(int i = N+1; i > 1; i--){
        v[i-1]= (s[i-1] - c[i-1]*v[i])/b[i-1];

    }
}

//specific gaussian elimination - sub and above diagonals have the same value = -1 -- cutting # of FLOPS with this "trick"
void specific_solver(int N, double* b, double* s, double* v){
    //forward
    for(int i = 1; i <= N+1; i++){
        b[i] = (i+1.0)/((double)i);
        s[i] += s[i-1]/(b[i-1]);

    }
    //backward
    v[N-1] = s[N-1]/b[N-1];

    for(int i = N+1; i > 1; i--){
        v[i-1] = (s[i-1] + v[i])/b[i-1];
    }
}

double rel_error(double sim, double ex){
    //assures division by 0
    double eps = 1e-6;
    if(ex <= eps){
        return 0.0;
    }
    else{
        return (abs((sim - ex)/(ex)));
    }
}
