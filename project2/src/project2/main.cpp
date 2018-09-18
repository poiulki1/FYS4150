#include <iostream>
#include <armadillo>
#include <iomanip>
#include <cmath>
#include <cstdlib>

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

    if(argc < 1){
        cout << "Please give a correct agrument" << endl;
        exit(1);
    }
    else{
        exponent = atoi(argv[1]);
    }

    test_max();
    test_jacobi();
    exit(1);

    //for(int i = 1; i <= exponent; i++){
    n = (int) pow(10.0,exponent); //loope for hver N, samme som proj. 1

    double start = 0.0;
    double stop = 5.0;
    double dt = (stop-start)/((float) n);


    double eps = pow(10.0, -8); //tolerance ~0

    arma::mat A = arma::zeros<arma::mat>(n,n);
    arma::mat R = arma::zeros<arma::mat>(n,n);

    double diag = 2.0/(dt*dt);
    double non_diag = -1.0/(dt*dt);

    unsigned int k, l;

    fill_tridiagonal_matrix(A, n, diag, non_diag);

    //start time here
    double max = 1; //bigger than eps

    int max_iter = 10000000;
    int iter = 0;

    while(max > eps && iter < max_iter){
        find_max_offdiag(A, k, l, n);
        max = fabs(A(l,k));
        jacobirotate(A, R, k, l, n);
        iter++;
    }

    cout << iter << " transformation(s) were needed to diagonalize matrix A" << endl;

    arma::mat exact_eig_values(n,n);

    analytic_eig(exact_eig_values, n, non_diag, diag);

    for(int bb = 0; bb < n; bb++){
        cout << "A: " << A(bb,bb) << "..." << "exact eig val: " << exact_eig_values(bb,bb) << endl;
    }
    return 0;
}
