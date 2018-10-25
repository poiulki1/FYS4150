#include <iostream>
using namespace std;
#include "solver.h"
#include "vec3.h"
#include "body.h"
#include <iomanip>
#include <fstream>


solver::solver(double ddt, quantities& sSystem):
    //arg(argument);
    dt(ddt)

{
    solarSystem = &sSystem;
}

void solver::forward_euler(quantities &system, int iter, double beta){
    bool condition;
    condition = false;

    system.calculate_system(condition, beta);
    system.calculate_potential_energy();
    for(int i = 0; i < system.list_of_bodies.size(); i++){
        body &current_body = system.list_of_bodies[i];

        current_body.position += current_body.velocity*dt;
        current_body.velocity += current_body.acceleration*dt;



    }
    condition = true;
    system.calculate_system(condition, beta);
    for(int i = 0; i < system.list_of_bodies.size(); i++){
        body &current_body = system.list_of_bodies[i];
        system.information[i].update(current_body, iter);
    }


}

void solver::velocity_verlet(quantities &system, int iter, double beta){
    bool condition;
    condition = false;

    system.calculate_system(condition, beta);
    system.calculate_potential_energy();
    for(int i = 0; i < system.list_of_bodies.size(); i++){
        body &current_body = system.list_of_bodies[i];
        current_body.velocity += 0.5*current_body.acceleration*dt;
        current_body.position += current_body.velocity*dt;

    }

    condition = true;
    system.calculate_system(condition, beta);
    for(int i = 0; i < system.list_of_bodies.size(); i++){
        body &current_body = system.list_of_bodies[i];
        current_body.velocity+= 0.5*current_body.acceleration*dt;
    }

    //if(iter%100 == 0){
        for(int i = 0; i < system.list_of_bodies.size(); i++){
            body &current_body = system.list_of_bodies[i];
            system.information[i].update(current_body, iter);
        }
    //}






}


