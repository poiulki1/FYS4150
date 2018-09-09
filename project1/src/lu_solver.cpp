#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <armadillo>
using namespace std;


inline double f(double x){return 100.0*exp(-10.0*x);}
inline double exact(double x){return 1.0 - (1 - exp(-10))*x - exp(-10*x);}
double rel_error(double, double);
void LU_decomposistion(int, arma::vec &, arma::vec &, arma::vec &, arma::vec &);


int main(int argc, char* argv[]){
    int N, algo, exponent;
    if(argc < 3){
        cout << "You need to give two arguments: \n first: N - number of mesh points \n second: algo: general - 0, spec - 1" << endl;
        exit(1);
    }
    else{
        exponent = atoi(argv[1]);
        algo = atoi(argv[2]);
    }

    ofstream file_error ("error.txt");
    ofstream file_solution ("solution.txt");

    for(int i = 1; i <= exponent; i++){
        N = (int) pow(10.0,i);
        file_solution << N << endl;

    double dt = 1.0/(N+1);

    arma::vec a = arma::zeros<arma::vec>(N+1);
    arma::vec b = arma::zeros<arma::vec>(N+2);
    arma::vec c = arma::zeros<arma::vec>(N+1);
    arma::vec v = arma::zeros<arma::vec>(N+2);




    for(int i = 0; i < N+1; i++){
        a(i) = -1.0;
        b(i) = 2.0;
        c(i) = -1.0;
    }

    b(N+1) = 2.0;
    LU_decomposistion(N,a,b,c,v);

    double* err = new double[N+1];
    for(int k = 1; k < N+2; k++){

        err[k] = rel_error(v[k], exact(dt*k));

     }
     double biggest_error = err[0];

     for(int j = 1; j < N; j++){
        if(abs(err[j]) >= abs(biggest_error)){
                //cout << err[j] << endl;
                biggest_error = err[j];
         }
      }

     file_error << setprecision(8) << left << setw(12) << biggest_error
                                     << left << setw(12) << N
                                     << endl;

     for(int idx = 0; idx <N+2; idx++){
         file_solution << setprecision(8) << left << setw(20) << v[idx]
                                         << left << setw(20) << dt*idx
                                         << endl;
     }
    }



}

void LU_decomposistion(int N, arma::vec& a, arma::vec& b, arma::vec& c, arma::vec& v){
    clock_t start_time, finish_time;

    double dt = 1.0/(N+1.0);
    arma::vec y = arma::zeros<arma::vec>(N);
    arma::vec x = arma::zeros<arma::vec>(N);
    arma::vec function = arma::zeros<arma::vec>(N);

    arma::mat A = arma::zeros<arma::mat>(N,N);
    arma::mat L, U;

    for(int i = 0; i < N; i++){
        function[i] = dt*dt*f(dt*i);
        A(i,i) = b(i);
        if(i+1 != N){
            A(i+1,i) = a(i);
            A(i,i+1) = c(i);


        }

    }


    start_time = clock();
    arma::lu(L,U,A);
    y = arma::solve(L,function,arma::solve_opts::no_approx);
    x = arma::solve(U,y,arma::solve_opts::no_approx);

    for(int j = 1; j < N+2; j++){
        v[j] = x[j-1];
    }
    v[0] = v[N+1] = 0;
    finish_time = clock();
    cout << "time used: "<< (double)(finish_time - start_time)/CLOCKS_PER_SEC << " sec " << " N = " << N << endl;

}



double rel_error(double sim, double ex){
    double eps = 1e-6;
    if(ex <= eps){
        return 0.0;
    }
    else{
        return abs((sim - ex)/(ex));
    }
}
