#ifndef EIGENVALUE_SOLVER_H
#define EIGENVALUE_SOLVER_H
#include <armadillo>


    void jacobi(arma::mat&, arma::vec&, int, int, double);
    void jacobirotate(arma::mat&, arma::mat&, unsigned int&, unsigned int&, int&);
    void find_max_offdiag(arma::mat&, unsigned int&, unsigned int&, int&);

    void fill_tridiagonal_matrix(arma::mat&, arma::vec&, int, double, double);
    void fill_tridiagonal_matrix_no_potential(arma::mat&, int, double, double);
    void initialize_tridiagonal_matrix(arma::mat&, int, double, double);
    void eigen_value_check_armadillo(arma::mat&);

    void analytic_eig(arma::vec&, int, double, double);


#endif // EIGENVALUE_SOLVER_H


