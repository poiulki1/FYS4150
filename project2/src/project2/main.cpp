#include <iostream>
#include <armadillo>
#include <iomanip>
#include <cmath>

//header files
#include "time.h"
#include "eigenvalue_solver.h"


using namespace std;


//void eigen_value_check_armadillo(arma::mat &);
void fill_tridiagonal_matrix(arma::mat&, int, double, double);
//void find_max_offdiag(arma::mat&, unsigned int&, unsigned int&, unsigned int&); //trenger k og l som argumenter i denne?
//void jacobirotate(arma::mat&, arma::mat&, unsigned int&, unsigned int&, unsigned int&);

int main(int argc, char* argv[]){

    int n, exponent;
    if(argc < 2){
        cout << "Please give a correct agrument" << endl;
        exit(1);
    }
    else{
        exponent = atoi(argv[1]);
    }

    //for(int i = 1; i <= exponent; i++){
    n = (int) pow(10.0,exponent); //loope for hver N, samme som proj. 1
    //double dt = (stop-start)/((float) n);

    double eps = pow(10.0, -9); //tolerance ~0

    arma::mat A = arma::zeros<arma::mat>(n,n);
    arma::mat R = arma::zeros<arma::mat>(n,n);
    double diag = 2.0;//(dt*dt);
    double non_diag = -1.0;//(dt*dt);
    int k, l;


    fill_tridiagonal_matrix(A, n, diag, non_diag);

    eigen_value_check_armadillo(A);
    //start time here
    int max = 0.01;
    int max_iter = 10;
    int iter = 1;
    while(max > eps && iter <= max_iter){
        cout << max << "..." << iter << endl;
        find_max_offdiag(A, k, l, n, max);
        jacobirotate(A, R, k, l, n);

        iter++;
    }
    cout << iter << "transformation were needed to diagonalize matrix A" << endl;

    test_jacobi();
    return 0;
}

void fill_tridiagonal_matrix(arma::mat& matrix, int n, double diag, double non_diag){
    // diag and non_diag, are the values we want to fill given matrix, which becomes tridiagonal symmetric matrix
    for(int i = 0; i < n; i++){
        matrix(i,i) = diag;
        if(i != n-1 ){
            matrix(i+1,i) = non_diag;
            matrix(i,i+1) = non_diag;
        }
    }
}
