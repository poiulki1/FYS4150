#include <iostream>
#include "solver.h"
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

int main()
{   int MC = 5000;
    int initial_S = 300;
    int initial_I = 100;
    int initial_R = 0;
    double a = 4;
    double b = 1;
    double c = 0.5;
    double d = 0.1;
    double d_i = 0.1;
    double e = 0.2;
    double f = 0;


    double final_time = 15;
    double start_time = 0;
    int N_population = 400;


    //cout << "Hello World!, this is an apocalypse" << endl;
    solver methods(MC,initial_S,initial_I,initial_R, a, b, c, d,d_i,e,f,N_population, final_time, start_time);

    methods.execute_solve(true, 5);

    //cout << "initial state " << " " << methods.S[0] << " end state " << " " << methods.S[MC] << endl;
    //cout << "initial state " << " " << methods.I[0] << " end state " << " " << methods.I[MC] << endl;
    //cout << "initial state " << " " << methods.R[0] << " end state " << " " << methods.R[MC] << endl;

    return 0;
}
