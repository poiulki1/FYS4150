#include "store_values.h"
#include <iostream>
#include <fstream>

store_values::store_values(int N)
{
    position.resize(int(N/100.0));
    //velocity.resize(N);
    //kinetic_energy.resize(N);
    //potential_energy.resize(N);
    //angular_momentum.resize(N);
    //momentum.resize(N);
    //theta.resize(N);



}

void store_values::update(body &noe, int i){

   position[i] = noe.position;
    //velocity[i] = noe.velocity;
    //kinetic_energy[i] = noe.K;
   //potential_energy[i] = noe.P;
    //angular_momentum[i] = noe.L;
    //momentum[i] = noe.moment;
    //theta[i] = noe.angle;
    double r;


}
void store_values::write_to_file(double beta){

    string mm = to_string(beta);

    ofstream file;
    //ofstream mfile;
    file.open(arg + ".txt");
    //mfile.open(arg + "_quantity_info.txt");


    for(int i= 0; i < position.size(); i++){

            //file << theta[i] << " ";
            file << position[i][0] << " " << position[i][1] << " " << position[i][2] << " ";
            //file << velocity[i][0] << " " << velocity[i][1] << " " << velocity[i][2] << " ";
            file << endl;

            //mfile << kinetic_energy[i] << " ";
            //mfile << potential_energy[i] << " ";
            //mfile << momentum[i] << " ";
            //mfile << angular_momentum[i] << " ";
            //mfile << endl;


    }

    file.close();
    //mfile.close();
    cout << "file written " << arg + ".txt" << endl;
    //cout << "file written " << arg + "_quantity_info.txt" << endl;
    cout << "_________________________________________________" << endl;
    double c = 0;
}
