//#include <iostream>
#include <armadillo>
#include "metropolis.h"

using namespace arma;
using namespace std;



int main(int argc, char *argv[])
{
    int n_spins, n_mc, n_temp;
    double temp_max, temp_min, E, M, comp_fac[17];
    ofstream output_file;
    vec quantities;
    mat s_matrix;

    // Handle arguments
    int arg_lim = 5;
    if(argc < arg_lim){
        cout << "Error: too few arguments" << endl;
        exit(1);
    }
    else{
        n_spins = stoi(argv[1]);
        n_mc = stoi(argv[2]);
        n_temp = stoi(argv[3]);
        temp_max = stod(argv[4]);
        temp_min = stod(argv[5]);
        output_file.open("data.txt");
    }

    output_file << setw(15) << "Temperature";
    output_file << setw(15) << "<E>";
    output_file << setw(15) << "<E2>";
    output_file << setw(15) << "var[E]";
    output_file << setw(15) << "Cv";
    output_file << setw(15) << "<M>";
    output_file << setw(15) << "<M2>";
    output_file << setw(15) << "var[M]";
    output_file << setw(15) << "Xi";
    output_file << endl;


    s_matrix = ones(n_spins, n_spins);

    double temp_step = (temp_max - temp_min)/n_temp;
    for(double temp = temp_min; temp <= temp_max; temp += temp_step){
        // Quantities vector will store expected values of energy and magnetization
        quantities = zeros(5);

        // Energy and Magnetization
        E = M = 0;

        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) comp_fac[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) comp_fac[de+8] = exp(-de/temp);

        // initializing for different temp
        initialize(n_spins, temp, s_matrix, E, M);

        for(int mc_step = 0; mc_step <= n_mc; mc_step++){
            // Metropolis algorithm
            metropolis(n_spins, s_matrix, E, M, comp_fac);
            // Update expectation values
            quantities[0] += E;
            quantities[1] += E*E;
            quantities[2] += M;
            quantities[3] += M*M;
            quantities[4] += fabs(M);
        }

        write_to_file(output_file, n_spins, n_mc, temp, quantities);
    }

    return 0;
}
