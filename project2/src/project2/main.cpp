//sjekk kommentarer i eigenvalue_solver.cpp under test_jacobi func

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <cmath>
#include <cstdlib>

//header files
#include "time.h"
#include "jacobi.h"
#include "test.h"
#include "quantum.h"

using namespace std;

int main(int argc, char* argv[]){

    int n, exponent, max_iter;

    if(argc < 2){
        cout << "Please give a correct agruments: power of 10 to determinate steps, max number of rotations" << endl;
        exit(1);
    }
    else{
        exponent = atoi(argv[1]);
        max_iter = atoi(argv[2]);

        if((exponent*max_iter) == 0){
            cout << "All of the arguments need to be numbers" << endl;
            exit(1);
        }
    }

    //unit tests
    test_max();
    test_jacobi();

    //for(int i = 1; i <= exponent; i++){
    n = (int) pow(10.0,exponent); //loope for hver N, samme som proj. 1

    double start = 0.0;
    double stop = 1.0;

    double eps = pow(10.0, -8); //tolerance ~0

    //memory allocation
    arma::mat A = arma::zeros<arma::mat>(n,n); //matrix to be initialized later
    arma::vec eig_val(n); //vector to store numerical eigenvalues - sorted A diagonal elements
    arma::vec exact_eig_values(n); //vector to store analytical eigenvalues

    //initialize_tridiagonal_matrix(A, n, start, stop);
    repulsive_hamilton(A, start, stop, n, 0.5);

    //analytical eigenvalues for tridiagonal matrix with the same elements on diagonal, and different but the same on off diag
    //analytic_eig(exact_eig_values, n, non_diag, diag);

    //jacobi algo. - rotating and finding the max off diag element repeated untill diagonalized or untill max iter rotations
    jacobi(A, eig_val, n, max_iter, eps);


    return 0;
}
