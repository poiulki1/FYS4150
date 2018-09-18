//sjekk kommentarer i eigenvalue_solver.cpp under test_jacobi func

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <cmath>
#include <cstdlib>

//header files
#include "time.h"
#include "eigenvalue_solver.h"


using namespace std;

int main(int argc, char* argv[]){

    int n, exponent;

    if(argc < 1){
        cout << "Please give a correct agrument" << endl;
        exit(1);
    }
    else{
        exponent = atoi(argv[1]);
    }

    test_max();
    test_jacobi();

    //for(int i = 1; i <= exponent; i++){
    n = (int) pow(10.0,exponent); //loope for hver N, samme som proj. 1

    double start = 0.0;
    double stop = 1.0;
    double dt = (stop-start)/((float) n);

    double eps = pow(10.0, -8); //tolerance ~0

    arma::mat A = arma::zeros<arma::mat>(n,n);
    arma::mat R = arma::eye<arma::mat>(n,n);

    arma::vec eig_val(n); //vector to store eigenvalues - sorted A diagonal elements

    double diag = 2.0/(dt*dt);
    double non_diag = -1.0/(dt*dt);

    unsigned int k, l;

    fill_tridiagonal_matrix(A, n, diag, non_diag);

    //start time here
    double max = 0.01; //bigger than eps

    int max_iter = 1000000;
    int iter = 0;

    while(max > eps && iter < max_iter){
        find_max_offdiag(A, k, l, n);
        max = fabs(A(l,k));
        jacobirotate(A, R, k, l, n);
        iter++;
    }

    cout << iter << " transformation(s) were needed to diagonalize matrix A" << endl;

    arma::vec exact_eig_values(n);

    analytic_eig(exact_eig_values, n, non_diag, diag);

    eig_val = arma::sort(diagvec(A));

    for(int bb = 0; bb < n; bb++){
        cout << "num: " << eig_val(bb) << "..." << "anal: " << exact_eig_values(bb) << endl;
    }
    return 0;
}
