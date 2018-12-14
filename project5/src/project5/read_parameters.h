#ifndef READ_PARAMETERS_H
#define READ_PARAMETERS_H
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace  std;

class read_parameters
{
public:
    read_parameters();

    static read_parameters* m_instance;

    static void parameters(std::string location);

    static int MC_cycles;
    static bool MC_cycles_set;

    static int dimension;
    static bool dimension_set;


    static int Npopulation;
    static bool Npopulation_set;

    static int Nsamples;
    static bool Nsamples_set;


    static int init_S;
    static bool init_S_set;

    static int init_I;
    static bool init_I_set;

    static int init_R;
    static bool init_R_set;

    static double a;
    static bool a_set;

    static double b;
    static bool b_set;

    static double c;
    static bool c_set;

    static double d;
    static bool d_set;

    static double d_i;
    static bool d_i_set;

    static double e;
    static bool e_set;

    static double f;
    static bool f_set;

    static double omega;
    static bool omega_set;

    static double max_diviation;
    static bool max_diviation_set;

    static double final_time;
    static bool final_time_set;

    static double start_time;
    static bool start_time_set;

    static string variance_save;
    static bool variance_save_set;

    static string filename_save;
    static bool filename_save_set;

    static bool solver;
    static bool solver_set;




    static bool numerical;
    static bool numerical_set;

};

#endif // READ_PARAMETERS_H


