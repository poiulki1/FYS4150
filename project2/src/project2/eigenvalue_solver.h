#ifndef EIGENVALUE_SOLVER_H
#define EIGENVALUE_SOLVER_H
#include <armadillo>



    void eigen_value_check_armadillo(arma::mat&);
    void analytic_eig(arma::mat&, int, double, double);
    void find_max_offdiag(arma::mat&, unsigned int&, unsigned int&, int&);
    void jacobirotate(arma::mat&, arma::mat&, unsigned int&, unsigned int&, int&);
    void fill_tridiagonal_matrix(arma::mat&, int, double, double);
    void test_jacobi();
    void test_max();


#endif // EIGENVALUE_SOLVER_H


