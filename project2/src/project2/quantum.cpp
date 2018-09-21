#include "quantum.h"
#include "jacobi.h"
#include <armadillo>

using namespace std;

void hamilton(arma::mat& A, double rho_min, double rho_max, int size){

    arma::vec potential(size);

    double h = (rho_max - rho_min)/((float) size);
    double non_diag = -1.0/(h*h);
    double diag = 2.0/(h*h);

    for(int i = 0; i < size; i++){
        potential(i) = pow(rho_min + i*h, 2);
    }

    fill_tridiagonal_matrix(A, potential, size, diag, non_diag);
}

void repulsive_hamilton(arma::mat& A, double rho_min, double rho_max, int size, double omega){

    arma::vec potential(size);

    double h = (rho_max - rho_min)/((float) size);
    double non_diag = -1.0/(h*h);
    double diag = 2.0/(h*h);

    double K; //constant
    potential(0) = h;
    for(int i = 1; i < size; i++){
        K = rho_min + i*h;
        potential(i) = (omega*omega*K*K) + (1.0/K);
    }

    fill_tridiagonal_matrix(A, potential, size, diag, non_diag);
}
