#include <iostream>
#include <math.h>
#include <armadillo>
using namespace std;

double acc(arma::vec&);
void Euler_solver(int, double, arma::mat&, arma::mat&, arma::vec&);

int main()
{

    return 1;
}

double acceleration(arma::vec& r){
    double pi = acos(-1);

    return ((4*pi*pi)/(arma::norm(r)));
}

void Euler_solver(int N, double dt, arma::mat& pos, arma::mat& vel, arma::vec& time){

    for(int i = 0; i<N; i++){
        pos(i+1) = pos(i) + vel(i)*dt;

        vel(i+1) = vel(i) + acceleration(pos(i))*dt;
        time(i+1) = time(i) + dt;
    }
}


//void Verlet_solver(int N, double dt){
//
//}
