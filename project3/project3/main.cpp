#include <iostream>
#include <math.h>
#include <armadillo>
#include "vec3.h"
#include "body.h"
#include "solver.h"
using namespace std;

//void write_to_file(char file_name, vec3 data);

int main()
{
    vec3 earth_pos(1.0, 0.0, 0.0);
    vec3 earth_vel(0.0, 3*M_PI, 0.0);

    vec3 sun_pos(0.20,0.0,0.0);
    vec3 sun_vel(0.0,0.0,0.0);

    body Earth(0.000003, earth_pos, earth_vel);
    body Sun(1.0, sun_pos, sun_vel);

    solver solar_system;

    solar_system.add_new_body(Sun);
    solar_system.add_new_body(Earth);

    int N = 10000;
    double final_time = 1.0;

    double dt = final_time/((double) N);


    cout << Earth.position << endl;
    for(int i=0; i<N; i++){

        solar_system.euler_one_step(dt);
    }

    cout << Earth.position << "finish" << endl;


    return 1;
}

/*
void write_to_file(){
    ofstream data (file_name);

}*/
