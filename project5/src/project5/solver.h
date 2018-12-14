#ifndef SOLVER_H
#define SOLVER_H

#include <functional>
#include <random>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include "savetofile.h"
using namespace std;
using namespace arma;

class solver
{
private:
    inline double dt1(double a_value, double N){
        return 4.0/((a_value)*N);
    }

    inline double dt2(double b_value, double N){
        return 1.0/(b_value*N);
    }

    inline double dt3(double c_value, double N){
        return 1.0/(c_value*N);
    }

    inline double dt4(double d_value, double N){
        return 1.0/(d_value*N);
    }

    inline double dt5(double dI_value, double d_value, double N){
        return 1.0/((dI_value)*N);
    }

    inline double dt6(double d_value, double N){
        return 1.0/(d_value*N);
    }

    inline double dt7(double e_value, double N){
        return 1.0/(e_value*N);
    }


    int ss, ii, rr, dd, dss, diii, drrr, bss;

    double S_to_I, I_to_R, R_to_S, S_to_D, I_to_Di,I_to_D, R_to_D, B_to_S, S_to_R;
    double S, I , R;
    double init_S, init_I, init_R;
    double susceptible, infected, recovered, born, death_sus, death_rec, death_infected;

    void update_a(double time);
    double a0, f0;
    void reset_matrix(mat &A, rowvec &B);
    void set_initial(bool Mc_arg, mat &M);
    void update_f(double time);
    void calculate_variance_and_sigma(int b, rowvec &time, mat &sampled_SIR, mat &unsampled_SIR,
                                      int nsamples, savetofile &save_obj);
public:


    solver(int Mc, double initial_S, double initial_I,double initial_R,
           int N_population, double final_time, double start_time);
    ~solver();

    int MCcycles, N;
    double ft,ts;
    int number_of_functions = 3;




    double tid;

    double *parameter;

    double s,i,r;
    double condition1, condition2, condition3;

    double k1,k2,k3,k4;



    void MonteCarlo(double dt);
    void execute_solve(string arg, bool Mc_arg, double B, int nsamples, savetofile &save_obj);

    void change_b(double value);

    rowvec derivative(int time, rowvec &y);
    void rk4(int time, mat &y, int dim, int index, double h);

    double dt(double N);

    double limit;
    double omega;


    void Npopulation(double S_, double I_, double R_);


    void initialize_parameters(double a, double b, double c, double d, double d_i, double e, double f, double deviation, double w);
};

#endif // SOLVER_H
