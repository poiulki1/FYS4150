#include <iostream>
#include "solver.h"
#include <vector>
#include <armadillo>
#include "read_parameters.h"
using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{

    read_parameters::parameters("../parameter/parameter.txt");


    solver run_solver(read_parameters::MC_cycles,read_parameters::init_S,read_parameters::init_I,
                      read_parameters::init_R, read_parameters::Npopulation,
                      read_parameters::final_time, read_parameters::start_time);

    savetofile save_obj("run", "a(t)_0.5_0_0_0_0.7");

    run_solver.initialize_parameters(read_parameters::a, read_parameters::b, read_parameters::c, read_parameters::d,
                                     read_parameters::d_i,read_parameters::e,read_parameters::f,
                                     read_parameters::max_diviation, read_parameters::omega);

    run_solver.execute_solve("mc",true , 5,read_parameters::Nsamples ,save_obj);
    run_solver.execute_solve("rk4",false , 5,read_parameters::Nsamples ,save_obj);

    return 0;
}
