#include "read_parameters.h"

read_parameters::read_parameters()
{

}


/*
Loads the parameters from a parameter txt file. This class is static, so that the parameters only need to be loaded
once. The parameters in then retrieved by using Parameters::name_of_float_vars.
In the txt file the parameters are given simply as
parameter_name parameter_value #Comment
'#' works as a marker for comments (a simple space between the value and comment may work
but use # to be safe.
The parameter name has to coinside with the names in this class. If unknow name is found in the
file, or a parameter if not set in the file, an error will be given.
*/

void read_parameters::parameters(string location){
    ifstream infile;
    string line;

    infile.open(location);
    if(infile.fail()){
        cout << "Could not find parameter file at " << location << endl;
        exit(EXIT_FAILURE);
    }
    infile.clear();
    infile.seekg(0,ios::beg);
    for (string l; getline(infile,l);)
    {
        stringstream ss(l);

        string name;
        double float_var;


        ss >> name >> float_var;

        if (name.front() == '#'){
            continue;
        }

        else if(name == "MC_cycles"){
            MC_cycles = float_var;
            if(MC_cycles>0){
                MC_cycles_set=true;
            }
        }

        else if(name == "Number_population"){
            Npopulation = float_var;
            if(Npopulation>0){
                Npopulation_set = true;
            }
        }

        else if(name == "Nsamples"){
            Nsamples = float_var;
            if(Nsamples>0){
                Nsamples_set = true;
            }
        }

        else if(name == "initial_susceptible"){
            init_S = float_var;
            if(init_S >= 0){
                init_S_set =true;
            }
        }


        else if(name == "initial_infected"){
            init_I = float_var;
            if(init_I >= 0){
                init_I_set =true;
            }
        }
        else if(name == "initial_recovered"){
            init_R = float_var;
            if(init_R >= 0){
                init_R_set =true;
            }
        }
        else if(name == "transmission_of_disease"){
            a = float_var;
            if(a >= 0){
                a_set =true;
            }
        }

        else if(name == "rate_of_recovery"){
            b = float_var;
            if(b >= 0){
                b_set =true;
            }
        }
        else if(name == "loss_of_immunity"){
            c = float_var;
            if(c >= 0){
                c_set =true;
            }
        }

        else if(name == "death_rate"){
            d = float_var;
            if(d >= 0){
                d_set =true;
            }
        }

        else if(name == "death_rate_infected"){
            d_i = float_var;
            if(d_i >= 0){
                d_i_set =true;
            }
        }

        else if(name == "birth_rate"){
            e = float_var;
            if(e >= 0){
                e_set =true;
            }
        }
        else if(name == "vaccine"){
            f = float_var;
            if(f >= 0){
                f_set =true;
            }
        }
        else if(name == "max_diviation"){
            max_diviation = float_var;
            if(max_diviation >= 0){
                max_diviation_set = true;
            }
        }

        else if(name == "omega"){
            omega = float_var;
            if(omega >= 0){
                omega_set = true;
            }
        }

        else if(name == "solver"){
            solver = (bool)float_var;
            solver_set = true;

        }


        else if(name == "final_time"){
            final_time = float_var;
            if(final_time > 0){
                final_time_set = true;
            }
        }

        else if(name == "start_time"){
            start_time = float_var;
            if(start_time >= 0){
                start_time_set = true;
            }
        }


        else{
            cout << "Unknown variable found: " << name << endl;
            exit(EXIT_FAILURE);
        }
    }

    if(!MC_cycles_set){
        cout << "MC cycles not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!init_S_set){
        cout << "Initial population for susceptible is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!init_I_set){
        cout << "Initial infected for susceptible is not set!" << endl;
        exit(EXIT_FAILURE);
    }
    else if(!init_R_set){
        cout << "Initial population for recovered is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!a_set){
        cout << "transmission rate of the infection is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!b_set){
        cout << "rate of recovery is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!c_set){
        cout << "loss of immunity is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!d_set){
        cout << "death rate is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!d_i_set){
        cout << "death rate of infected is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!e_set){
        cout << "birth rate is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!f_set){
        cout << "vaccine rate is not set!" << endl;
        exit(EXIT_FAILURE);
    }
    else if(!max_diviation_set){
        cout << "max diviation from constant a is not set!" << endl;
        exit(EXIT_FAILURE);
    }
    else if(!omega_set){
        cout << "omega is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!solver_set){
        cout << "Type of solver is not given!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!final_time_set){
        cout << "final time is not set!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!start_time_set){
        cout << "start time is not given!" << endl;
        exit(EXIT_FAILURE);
    }

    else if(!Nsamples_set){
        cout << "Number of sampling is not defined!" << endl;
        exit(EXIT_FAILURE);
    }


}



bool read_parameters::MC_cycles_set=false;
bool read_parameters::Npopulation_set = false;
bool read_parameters::init_S_set=false;
bool read_parameters::init_I_set=false;
bool read_parameters::init_R_set=false;
bool read_parameters::a_set=false;
bool read_parameters::b_set=false;
bool read_parameters::c_set=false;
bool read_parameters::d_set=false;
bool read_parameters::d_i_set=false;
bool read_parameters::e_set=false;
bool read_parameters::f_set=false;
bool read_parameters::max_diviation_set=false;
bool read_parameters::omega_set=false;
bool read_parameters::final_time_set=false;
bool read_parameters::start_time_set=false;
bool read_parameters::solver_set=false;
bool read_parameters::Nsamples_set=false;

int read_parameters::MC_cycles=0;
int read_parameters::Npopulation = 0;
int read_parameters::Nsamples = 0;
int read_parameters::init_S=0;
int read_parameters::init_I=0;
int read_parameters::init_R=0;
double read_parameters::a=0;
double read_parameters::b=0;
double read_parameters::c=0;
double read_parameters::d=0;
double read_parameters::d_i=0;
double read_parameters::e=0;
double read_parameters::f=0;
double read_parameters::max_diviation=0;
double read_parameters::final_time=0;
double read_parameters::start_time=0;
double read_parameters::omega=0;
bool read_parameters::solver=false;

