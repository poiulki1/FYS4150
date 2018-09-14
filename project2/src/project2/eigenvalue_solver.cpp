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

void analytic_eig(arma::vec& eig_values, int dim, double a, double d){

    for(int i = 0; i < dim; i++){
        eig_values[i] = d + 2*a*cos(((i+1)*acos(-1))/(dim + 1.0));
    }
}




void find_max_offdiag(arma::mat& A, unsigned int& k, unsigned int& l, unsigned int& n, double& max){
  // Finds maximum element on upper diagonal of a matrix A
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = i+1; j < n; j++){
      if (fabs(A(i,j)) > max){
        max = fabs(A(i,j));
        l = i;
        k = j;
      }
    }
  }
  //return max;
}



void jacobirotate(arma::mat& A, arma::mat& R, unsigned int& k, unsigned int& l, unsigned int& n){
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

void test_jacobi(){
    unsigned int n = 3;

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
    int iter, max_iter;
    iter = 1;  max_iter = 3;
    double max = 0.01;

    double eps = pow(10.0, -8);

    while(max > eps && iter <= max_iter){
        find_max_offdiag(A,k,l,n, max);
        jacobirotate(A, R, k, l, n);
        iter++;
    }

    if(abs(A(0,0) - eig_val(0)) < eps && abs(A(1,1) - eig_val(1))< eps && abs(A(2,2) - eig_val(2))< eps){
        cout << "Jacobi test passed" << endl;
    }
    else{
        cout << "fail"<< endl;
        exit(1);
    }
}
