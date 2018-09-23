#include <armadillo>
#include "test.h"
#include "jacobi.h"

using namespace std;

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
        //PASSED
    }
    else{
        cout << "Finding max offdiag element test FAILED - need to check the find_max_offdiag function" << endl;
        exit(1);
    }

}

void test_jacobi(){
    int n = 4;
    unsigned long long size = 4;
    int index1, index2;

    arma::mat A = arma::zeros<arma::mat>(size,size);
    arma::mat R = arma::eye<arma::mat>(size,size);

    arma::vec exact_eig_val(size);
    arma::vec eig_val(size);

    A(0,0) = A(1,1) = A(2,2) = A(3,3) = 3;

    A(0,1) = A(1,0) = 1;
    A(1,2) = A(2,1) = 4;
    A(2,3) = A(3,2) = 5;

    exact_eig_val(0) = -3.4340;
    exact_eig_val(1) = 2.2229;
    exact_eig_val(2) = 3.7771;
    exact_eig_val(3) = 9.4340;

    unsigned int k, l;
    int iter, max_iter;
    iter = 1;  max_iter = 10000;


    double eps = pow(10.0, -8);
    double max = 0.1;

    while(max > eps && iter < max_iter){
        find_max_offdiag(A,k,l,n);
        max = fabs(A(l,k));
        jacobirotate(A, R, k, l, n);
        iter++;
    }

    eig_val = arma::sort(diagvec(A)); //sort eigenvalues and store in a vector


    //diagonalization/eigenvalue test
    if((fabs(eig_val(0) - exact_eig_val(0))) && (fabs(eig_val(1) - exact_eig_val(1))) && (fabs(eig_val(2) - exact_eig_val(2))) && (fabs(eig_val(3) - exact_eig_val(3)))){
        //PASSED
    }
    else{
        cout << "Jacobi diagonalization test FAILED - need to fix jacobi rotate function"<< endl;
        exit(1);
    }

    //orthogonality of eigenvectors test
    arma::vec vec1(size);
    arma::vec vec2(size);

    for(int ind = 0; ind < n; ind++){
        vec1(ind) = R(ind,0);
        vec2(ind) = R(ind,1);
    }

    if((arma::dot(vec1,vec2)) < eps){
        //PASSED
    }
    else{
        cout << "Jacobi orthogonality test FAILED - need to fix jacobi rotate function" << endl;
        exit(1);
    }

    //symmetry test
    index1 = rand() % n;
    index2 = rand() % index1-1;

    if(fabs(A(index1, index2) - A(index2, index1)) < eps){
        //PASSED
    }
    else{
        cout << "Symmetry of a matrix is gone after the Jacobi rotation(s)" << endl;
        exit(1);
    }
}
