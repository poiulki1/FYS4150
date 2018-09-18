#include "eigenvalue_solver.h"

using namespace std;

void eigen_value_check_armadillo(arma::mat& matrix){
        arma::vec eigenvalues;
        arma::mat eigenvectors;

        arma::eig_sym(eigenvalues, eigenvectors, matrix); //hvorfor ikke eig_sym?
        //arma::eig_sys(eigenvalues, matrix);
        //return eigenvalues

        //cout << eigenvalues << endl;
}

void analytic_eig(arma::mat& eig_values, int dim, double a, double d){

    for(int i = 0; i < dim; i++){
        eig_values(i,i) = d + 2*a*cos(((i+1)*acos(-1))/(dim + 1.0));
        //cout << eig_values(i) << endl;
    }
}


void find_max_offdiag(arma::mat& A, unsigned int& k, unsigned int& l, int& n){
  // Finds maximum element on upper diagonal of a matrix A
  double max = 0.0;
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = i+1; j < n; j++){
      if (fabs(A(i,j)) > max){
        max = fabs(A(i,j));
        l = i;
        k = j;
      }
    }
  }
}


void jacobirotate(arma::mat& A, arma::mat& R, unsigned int& k, unsigned int& l, int& n){
  // Rotate: Find values of cos and sin
  double s,c;
  if(A(k,l) != 0.0){
    double t,tau;
    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    if(tau>0){
      t = 1.0/(tau + sqrt(1.0 + tau*tau)); //instead of -tau ± sqrt(1+tau*tau), make sure num. precision won't couse error
    }
    else{
      t = -1.0/(-tau + sqrt(1.0 + tau*tau));
    }

    c = 1/sqrt(1+t*t);
    s = c*t;
  }
  else{
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  // Changing the matrix elements with indices k and l
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0;
  A(l,k) = 0.0;
  // Then we change the remaining elements
  for(unsigned int i = 0; i<n; i++){
    if( i != k && i != l){
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
  r_ik = R(i,k);
  r_il = R(i,l);
  R(i,k) = c*r_ik - s*r_il;
  R(i,l) = c*r_il + s*r_ik;
  }
}


void test_max(){
    int n = 4;
    unsigned long long size = 4;

    double eps = pow(10,-8);

    arma::mat A = arma::zeros<arma::mat>(size,size);

    unsigned int k, l;
    double biggest_element = 1240.2;

    A(0,0) = A(1,1) = A(2,2) = A(3,3) = 3;
    A(2,2) = 9;
    A(0,1) = A(1,0) = 1;
    A(1,2) = A(2,1) = biggest_element;
    A(2,3) = A(3,2) = 9;

    find_max_offdiag(A,k,l,n);


    if((abs((A(l,k))-(abs(biggest_element))) < eps) && (k == 2) && (l == 1)){
        cout << "Finding max offdiag element test PASSED" << endl;
    }
    else{
        cout << "Finding max offdiag element test FAILED - need to check the find_max_offdiag function" << endl;
        exit(1);
    }

}

void test_jacobi(){
    int n = 4;
    unsigned long long size = 4;
    arma::mat A = arma::zeros<arma::mat>(size,size);
    arma::mat R = arma::zeros<arma::mat>(size,size);

    arma::vec eig_val(size);

    A(0,0) = A(1,1) = A(2,2) = A(3,3) = 3;

    A(0,1) = A(1,0) = 1;
    A(1,2) = A(2,1) = 4;
    A(2,3) = A(3,2) = 5;

    eig_val(0) = -3.4340;
    eig_val(1) = 2.2229;
    eig_val(2) = 3.7771;
    eig_val(3) = 9.4340;

    unsigned int k, l;
    int iter, max_iter;
    iter = 1;  max_iter = 10000;


    double eps = pow(10.0, -8);
    double max = 0.1;
    A.print();
    while(max > eps && iter < max_iter){
        find_max_offdiag(A,k,l,n);
        max = fabs(A(l,k));
        jacobirotate(A, R, k, l, n);
        iter++;
    }
    A.print();
    exit(1);

    if((abs(A(0,0) - eig_val(0)) < eps) && (abs(A(1,1) - eig_val(1))< eps) && (abs(A(2,2) - eig_val(2))< eps) && (abs(A(3,3) - eig_val(3))< eps)){
        cout << "Jacobi test PASSED" << endl;
    }
    else{
        cout << "Jacobi test FAILED - need to check jacobi rotate function"<< endl;
        exit(1);
    }
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
