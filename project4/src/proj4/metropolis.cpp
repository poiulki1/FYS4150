#include "metropolis.h"


void initialize(int n_spins, double temp, mat& s_matrix, double& E, double& M){
    for (int x_ind = 0; x_ind < n_spins; x_ind++){
        for (int y_ind = 0; y_ind < n_spins; y_ind++){
            M += (double) s_matrix(x_ind,y_ind);
            E -= (double) (s_matrix(x_ind, y_ind)*(s_matrix(periodic(x_ind-1, n_spins),y_ind) + s_matrix(x_ind,periodic(y_ind-1, n_spins))));
        }
    }
}


void metropolis(int n_spins, mat& s_matrix, double& E, double& M, double *comp_fac){

    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    for(int i = 0; i < n_spins; i++){
        for(int j = 0; j < n_spins; j++){
            int x = int(dist(gen)*double(n_spins));
            int y = int(dist(gen)*double(n_spins));

            int dE = 2*s_matrix(x,y)*(s_matrix(periodic(x+1,n_spins), y) + s_matrix(periodic(x-1,n_spins), y) + s_matrix(x, periodic(y+1,n_spins)) + s_matrix(x, periodic(y-1,n_spins)));

            if(dist(gen) <= comp_fac[int(dE)+8]){
                s_matrix(x,y) *= -1;
                E += (double) dE;
                M += (double) 2*s_matrix(x,y);
            }
        }
    }
}


void write_to_file(ofstream &file_name, int n_spins, int n_mc, double temperature, vec arr){
    double mc_norm = 1/((double) (n_mc)); // divided by total number of cycles
    double s_norm = 1.0/(n_spins*n_spins);

    double exp_E = arr[0]*mc_norm;
    double exp_E2 = arr[1]*mc_norm;
    double exp_M = arr[2]*mc_norm;
    double exp_M2 = arr[3]*mc_norm;
    double exp_Mabs = arr[4]*mc_norm;
    // all expectation values are per spin, multiply by s_norm = 1/(n_spins*n_spins)
    double var_E = (exp_E2 - exp_E*exp_E)*s_norm;
    double var_M = (exp_M2 - exp_M*exp_M)*s_norm;

    double var_Mabs = (exp_M2 - exp_Mabs*exp_Mabs)*s_norm; // Using this only in xi calculation

    double Cv = var_E/(temperature*temperature);
    double xi = var_Mabs/temperature;

    file_name << setiosflags(ios::showpoint | ios::uppercase);
    file_name << setw(15) << setprecision(8) << temperature;
    file_name << setw(15) << setprecision(8) << exp_E*s_norm;
    file_name << setw(15) << setprecision(8) << exp_E2*s_norm;
    file_name << setw(15) << setprecision(8) << var_E;
    file_name << setw(15) << setprecision(8) << Cv;
    file_name << setw(15) << setprecision(8) << exp_M*s_norm;
    file_name << setw(15) << setprecision(8) << exp_M2*s_norm;
    file_name << setw(15) << setprecision(8) << var_M;
    file_name << setw(15) << setprecision(8) << xi;
    file_name << endl;
} // end output function
