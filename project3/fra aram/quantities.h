#ifndef QUANTITIES_H
#define QUANTITIES_H
#include <cmath>
#include <vector>
#include <vec3.h>
#include "body.h"
#include "store_values.h"
#include <fstream>
using namespace std;

class quantities
{
public:
    quantities();

    friend class body;

    //vectors
    vector<body> list_of_bodies;
    vector<store_values> information;




    //vec3 variables
    vec3 a;
    vec3 R_CoM;




    //constants
    double g = 4*M_PI*M_PI;
    double c = 63239.7263;

    //variable
    int total_bodies;
    double total_kinetic_energy;
    double total_potential_energy;
    double l;
    double distance;
    double F;







    //functions
    void add_new_body(body new_body, store_values store_info);
    void calculate_system(bool condition, double beta);
    double force(body current_body, body other_body, double R);
    double kinetic_energy(body current_body);
    double beta_force(body current_body,body other_body, double distance, double beta);
    double potential_energy(body current_body, body other_body, double distance);
    double single_momentum(body current_body);
    double angular_momentum(body current_body);
    void write_to_file(double beta);
    void COM();
    void calculate_potential_energy();

    //double total_energy(double K, double P);

    double relativistic_force(body current_planet, body other_planet, double distance);

};

#endif // QUANTITIES_H
