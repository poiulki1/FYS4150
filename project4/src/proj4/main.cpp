//#include <iostream>
#include <armadillo>
#include <mpi.h>
#include "metropolis.h"

using namespace arma;
using namespace std;


int main(int argc, char *argv[])
{
    int n_spins, n_mc, n_temp, my_rank, numprocs;
    double temp_max, temp_min, E, M, comp_fac[17];
    ofstream output_file, output_file_expE_temp1, output_file_expE_temp24;
    vec quantities_local, quantities_global;
    mat s_matrix;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Handle arguments
    int arg_lim = 6;

    if(my_rank == 0 && argc < arg_lim){
        cout << "Error: too few arguments! \n";
        cout << "Give: number of spins, number of mc cycles, number of temperature steps, max & min temperature";
        cout << "and output file name" << endl;
        exit(1);
    }
    if(my_rank == 0){
        // Initializing values needed to calculations, given from cml arguments
        n_spins = stoi(argv[1]);
        n_mc = stoi(argv[2]);
        n_temp = stoi(argv[3]);
        temp_max = stod(argv[4]);
        temp_min = stod(argv[5]);
        output_file.open(argv[6]);
        output_file_expE_temp1.open("4C_temp1_random.txt");
        output_file_expE_temp24.open("4C_temp24_random.txt");
    }

    if(my_rank == 0){
        // Writing the first line in data file (only for orientation while developing)
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
    }

    // Brodcast shared-variables to all the nodes
    MPI_Bcast(&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_mc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Partitioning the MC-cycles among the nodes
    int n_partitions = n_mc/numprocs;
    int mc_local_start = (my_rank * n_partitions) + 1;
    int mc_local_stop = (my_rank+1) * n_partitions;
    if(my_rank == numprocs - 1 && mc_local_stop < n_mc) mc_local_stop = n_mc;


    // Initialize the spin matrix
    //init_matrix(s_matrix, n_spins, "not_random");
    init_matrix(s_matrix, n_spins, "random");


    // Finding the temperature step and entering the temperature loop
    // Inside the loop, initialize the E and M, additionaly precalculating comp_fac = delta E needed in the Metropolis algo
    double temp_step = (temp_max - temp_min)/n_temp;
    for(double temp = temp_min; temp <= temp_max; temp += temp_step){
        // Quantities vector will store expected values of energy and magnetization
        quantities_global = zeros(5);
        quantities_local = zeros(5);

        // Energy and Magnetization
        E = M = 0;

        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) comp_fac[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) comp_fac[de+8] = exp(-de/temp);

        // initializing for different temp
        initialize(n_spins, temp, s_matrix, E, M);

        for(int mc_step = 0; mc_step <= n_mc; mc_step++){ // Without MPI
        //for(int mc_step = mc_local_start; mc_step < mc_local_stop; mc_step++){
            // Metropolis algorithm
            metropolis(n_spins, s_matrix, E, M, comp_fac);
            // Update expectation values
            quantities_local[0] += E;
            quantities_local[1] += E*E;
            quantities_local[2] += M;
            quantities_local[3] += M*M;
            quantities_local[4] += fabs(M);
            // Writing to file each 100 step, first column=mc steps, second=expected E, third=expected |M| (task 4c)
            if(mc_step != 0 && mc_step % 100 == 0 && temp == 1.0){
                output_file_expE_temp1 << setw(15) << mc_step;
                output_file_expE_temp1 << setw(15) << setprecision(8) << (quantities_local[0]/((double) mc_step*n_spins*n_spins));
                output_file_expE_temp1 << setw(15) << setprecision(8) << (quantities_local[4]/((double) mc_step*n_spins*n_spins));
                output_file_expE_temp1 << endl;
            }
            if(mc_step != 0 && mc_step % 100 == 0 && temp == 2.4){
                output_file_expE_temp24 << setw(15) << mc_step;
                output_file_expE_temp24 << setw(15) << setprecision(8) << (quantities_local[0]/((double) mc_step*n_spins*n_spins));
                output_file_expE_temp24 << setw(15) << setprecision(8) << (quantities_local[4]/((double) mc_step*n_spins*n_spins));
                output_file_expE_temp24 << endl;
            }
        }
        /*
        // quantities are armadillo vectors/double pointer, index [0] tells about the start point only
        MPI_Reduce(&quantities_local[0], &quantities_global[0], 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        //cout << "local: " << quantities_local << "from rank: " << my_rank << "temp: " << temp << endl;

        // Write expectation values to file for the according temperature
        if(my_rank == 0){
            //cout << "hei fra procs nr " << my_rank << endl;
            //cout << "global: " << quantities_global << "temp: " << temp << endl;
            write_to_file(output_file, n_spins, n_mc, temp, quantities_global);
        }*/
    }
    MPI_Finalize();
    return 0;
}
