#ifndef EIGENVALUE_SOLVER_H
#define EIGENVALUE_SOLVER_H
#include <armadillo>



    void eigen_value_check_armadillo(arma::mat&);
    void analytic_eig(arma::vec&, int, double, double);
    void find_max_offdiag(arma::mat&, unsigned int&, unsigned int&, unsigned int&, double&);
    void jacobirotate(arma::mat&, arma::mat&, unsigned int&, unsigned int&, unsigned int&);
    void test_jacobi();


#endif // EIGENVALUE_SOLVER_H


