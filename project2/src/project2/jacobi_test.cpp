#include <iostream>
#include <armadillo>
#include <cmath>

#include "eigenvalue_solver.h"

/*
int main(){

    int n = 3;

    arma::mat A(n,n);
    arma::mat R(n,n);

    arma::vec eig_val(n);

    A(0,0) = A(1,1) = 3;
    A(2,2) = 5;
    A(0,1) = A(1,0) = 1;
    A(0,2) = A(1,2) = A(2,0) = A(2,1) = -1;

    eig_val(0) = 2.0;
    eig_val(1) = 3.0;
    eig_val(2) = 6.0;

    unsigned int k, l;
    int iter, max, max_iter;
    iter = 1; max = 3; max_iter = 10;

    double eps = pow(10.0, -9);

    while(max > eps && iter <= max_iter){
        find_max_offdiag(A,k,l,n);
        jacobirotate(A, R, k, l, n);
        iter++;
    }

    if(abs(A(0,0) - eig_val(0)) < eps && abs(A(1,1) - eig_val(1))< eps && abs(A(2,2) - eig_val(2))< eps){
        cout << "Jacobi test passed" << endl;
    }
    else{
        exit(1);
    }
    return 0;
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
    double dt = (stop-start)/((float) n);

    double eps = pow(10.0, -9); //tolerance ~0

    arma::mat A = arma::zeros<arma::mat>(n,n);
    double diag = 2.0/(dt*dt);
    double non_diag = -1.0/(dt*dt);
    int k, l, iter;


    fill_tridiagonal_matrix(A, n, diag, non_diag);

    eigen_value_check_armadillo(A);
    //start time here

    while(max > eps && iter <= max_iter){
        max = 5.0;

        find_max_offdiag(A, k, l, n); //trenger vi k og l som arg? - denne funk oppdaterer max og k og l
        jacobirotate(A, R, k, l, n);

        iter++;
    }
    cout << iter << "transformation were needed to diagonalize matrix A" << endl;

}*/
