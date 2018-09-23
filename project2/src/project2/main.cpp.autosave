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
    //handling arguments
    int n, matrix, algo;

    if(argc < 3){
        cout << "Please give a correct agruments" << endl;
        cout << "1) number of mesh-points" << endl;
        cout << "2) conditions: 1=without potential, 2=with harmonic potential (1 electron), 3=with harmonic potential and repulsiv forces (2 electrons)" << endl;
        cout << "3) algo: 1=jacobi, 2=armadillo" << endl;
        exit(1);
    }
    else{
        n = atoi(argv[1]);
        matrix = atoi(argv[2]);
        algo = atoi(argv[3]);

        if((n*matrix*algo) == 0){
            cout << "All of the arguments need to be numbers" << endl;
            exit(1);
        }
    }

    //unit tests
    cout << "Unit tests:" << endl;
    test_max();
    test_jacobi();

    cout << "All of the test PASSED" << endl;
    cout << "." << endl;
    cout << "." << endl;

    //rho_min and rho_max
    double start = 0.0;
    double stop = 5.0; //change manually

    double eps = pow(10.0, -8); //tolerance ~0
    int max_iter = 10e6; //maximum # of rotations/transformations

    //memory allocation
    arma::mat A = arma::zeros<arma::mat>(n,n); //matrix to be initialized later
    arma::mat R = arma::eye<arma::mat>(n,n); //matrix to store eigenvectors

    arma::vec eig_val(n); //vector to store numerical eigenvalues - sorted A diagonal elements

    arma::vec exact_eig_values(n); //vector to store analytical eigenvalues

    //declearing some time variables to calculate time used by different algorithms
    clock_t start_time, finish_time;
    double time_used;

    if(algo == 1){ //Jacobi
        if(matrix == 1){
            cout << "Initializing matrix A without any potential and solving using Jacobi transformation" << endl;
            initialize_tridiagonal_matrix(A, n, start, stop);
            //jacobi algo. - rotating and finding the max off diag element repeated untill diagonalized or untill max iter rotations
            start_time = clock();
            jacobi(A, R, eig_val, n, max_iter, eps);
            finish_time = clock();
            time_used = ((double) (finish_time - start_time))/CLOCKS_PER_SEC;
            cout << "Finding eigenpars with Jacobi transformation took: " << time_used << "seconds. With N=" << n << endl;
        }
        else if(matrix == 2){
            cout << "Initializing matrix A with a harmonic oscillator potential and solving using Jacobi transformation" << endl;
            hamilton(A, start, stop, n);
            //jacobi algo. - rotating and finding the max off diag element repeated untill diagonalized or untill max iter rotations
            jacobi(A, R, eig_val, n, max_iter, eps);
        }
        else{
            cout << "Initializing matrix A with a harmonic oscillator potential and Coulomb repulsiv forces, then solving using Jacobi transformation" << endl;
            double omega = 1; //change manually
            cout << "Omega: " << omega << endl;
            repulsive_hamilton(A, start, stop, n, omega);
            //jacobi algo. - rotating and finding the max off diag element repeated untill diagonalized or untill max iter rotations
            jacobi(A, R, eig_val, n, max_iter, eps);
        }
    }
    else{ //armadillo
        if(matrix == 1){
            cout << "Initializing matrix A without any potential and solving using armadillo function" << endl;
            initialize_tridiagonal_matrix(A, n, start, stop);
            start_time = clock();
            eigensolver_armadillo(A, eig_val, R);
            finish_time = clock();
            time_used = ((double) (finish_time - start_time))/CLOCKS_PER_SEC;
            cout << "Finding eigenpars with armadillo took: " << time_used << "seconds. With N=" << n << endl;

        }
        else if(matrix == 2){
            cout << "Initializing matrix A with a harmonic oscillator potential, then solving using armadillo function" << endl;
            hamilton(A, start, stop, n);
            eigensolver_armadillo(A, eig_val, R);
        }
        else{
            cout << "Initializing matrix A with a harmonic oscillator potential and Coulomb repulsiv forces, then solving using armadillo function" << endl;
            double omega = 0.01; //change manually
            cout << "Omega = " << omega << endl;
            repulsive_hamilton(A, start, stop, n, omega);
            eigensolver_armadillo(A, eig_val, R);
        }
    }



    /*
    cout << "eigenvalues" << endl;
    for(int kk = 0; kk < 5; kk++){
        cout << eig_val(kk) << endl;
    }*/
    //saving to binary file - changing the name manually
    eig_val.save("eig_val_rep_w001.bin", arma::raw_binary);
    R.save("eig_vec_rep_w001_max5.bin", arma::raw_binary);

    return 0;
}
