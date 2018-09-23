#include "jacobi.h"
#include "time.h"
#include <armadillo>

using namespace std;

void jacobi(arma::mat& A, arma::mat& R, arma::vec& eigenvalues, int size, int max_iter, double eps){
    //This function is applying jacobi rotation together with finding max offdiag element to diagonalize matrix A
    //In addition we are calculating the time used for the diagonalization

    unsigned int k, l;

    int iter = 0;
    double max = eps*100; //greater than eps

    while(max > eps && iter < max_iter){
        find_max_offdiag(A, k, l, size);
        max = fabs(A(l,k));
        jacobirotate(A, R, k, l, size);
        iter++;
    }

    eigenvalues = arma::sort(diagvec(A));

    cout << iter << " tranfsormation(s) was/were needed to diagonalize given matrix" << endl;

}

void find_max_offdiag(arma::mat& A, unsigned int& k, unsigned int& l, int& n){
  //This function finds maximum element on upper diagonal of a matrix A
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
  //This function is rotating the matrix A using jacobi method

  //Rotate: Find values of cos and sin
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
  //Changing the matrix elements with indices k and l
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0;
  A(l,k) = 0.0;
  //Then we change the remaining elements
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

void initialize_tridiagonal_matrix(arma::mat& A, int size, double start, double stop){
    double h = (stop-start)/((float) size);

    double diag = 2.0/(h*h);
    double non_diag = -1.0/(h*h);

    fill_tridiagonal_matrix_no_potential(A, size, diag, non_diag);
}

void fill_tridiagonal_matrix_no_potential(arma::mat& matrix, int n, double diag, double non_diag){
    // diag and non_diag, are the values we want to fill given matrix, which becomes tridiagonal symmetric matrix
    for(int i = 0; i < n; i++){
        matrix(i,i) = diag;
        if(i != n-1 ){
            matrix(i+1,i) = non_diag;
            matrix(i,i+1) = non_diag;
        }
    }
}

void fill_tridiagonal_matrix(arma::mat& matrix, arma::vec& potential, int n, double diag, double non_diag){
    // diag and non_diag, are the values we want to fill given matrix, which becomes tridiagonal symmetric matrix
    for(int i = 0; i < n; i++){
        matrix(i,i) = diag + potential(i);
        if(i != n-1 ){
            matrix(i+1,i) = non_diag;
            matrix(i,i+1) = non_diag;
        }
    }
}

void analytic_eig(arma::vec& eig_values, int dim, double a, double d){
    //This function fills eig_values vector with analytical eigenvalues for tridiagonal matrix

    for(int i = 0; i < dim; i++){
        eig_values(i) = d + 2*a*cos(((i+1)*acos(-1))/(dim + 1.0));
        //cout << eig_values(i) << endl;
    }
}


void eigensolver_armadillo(arma::mat& matrix, arma::vec& eigenvalues, arma::mat& eigenvectors){
        arma::eig_sym(eigenvalues, eigenvectors, matrix);

}
