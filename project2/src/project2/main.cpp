#include <iostream>
#include <armadillo>

using namespace std;

void eigen_value_check_armadillo(arma::mat &);
int main(int argc, char* argv[]){
    int n;
    if(argc < 2){
        cout << "Please give a correct agrument" << endl;
        exit(1);
    }
    else{
        n = atoi(argv[1]);
    }
    arma::mat A = arma::zeros<arma::mat>(n,n);
    //h^2 skal henge på matrise elementene denne gangen
    for(int i = 0; i < n; i++){
        A(i,i) = 2;
        if(i != n-1 ){
            A(i+1,i) = -1;
            A(i,i+1) = -1;

        }

    }

    eigen_value_check_armadillo(A);

}
void eigen_value_check_armadillo(arma::mat& matrix){
    arma::cx_vec eigenvalues;
    arma::cx_mat eigenvectors;

    arma::eig_gen(eigenvalues, eigenvectors, matrix);
    cout << eigenvalues << endl;
}


void jacobirotate(arma::mat& A, arma::mat& R, unsigned int& k, unsigned int& l, unsigned int& n){
  // Rotate: Find values of cos and sin
  double s,c;
  if(A(k,l) != 0.0){
    double t,tau;
    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    if(tau>0){
      t = 1.0/(tau + sqrt(1.0 + tau*tau)); //riktig? s.217 - skal man dele?
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

void find_max_offdiag(arma::mat& A, unsigned int& k, unsigned int& l, unsigned int& n){
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



//hei2
