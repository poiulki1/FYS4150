#include "quantities.h"
#include <iostream>

quantities::quantities()
{
}

void quantities::add_new_body(body new_body, store_values store_info)
{
    list_of_bodies.push_back(new_body);
    store_info.arg = new_body.argument;
    information.push_back(store_info);
}

double quantities::force(body current_body,body other_body, double R){

        F =  -(g*other_body.mass*current_body.mass)/(R*R);

    return F;
}

void quantities::calculate_system(bool condition, double beta){

    for(int i = 0; i < list_of_bodies.size(); i++){
        body &current_planet = list_of_bodies[i];
        a.zeros();

        for(int j = 0; j < list_of_bodies.size(); j++){
            if(i == j){
                continue;
            }

            body &other_planet = list_of_bodies[j];
            vec3 relative_coordinates = current_planet.position - other_planet.position;
            vec3 unit_vector = relative_coordinates/relative_coordinates.length();
            distance = relative_coordinates.length();



            a += force(current_planet, other_planet, distance)*unit_vector/current_planet.mass;
            //a += beta_force(current_planet, other_planet,distance, beta)*unit_vector/current_planet.mass;
            //a += relativistic_force(current_planet, other_planet, distance)*unit_vector/current_planet.mass;



        }

        current_planet.acceleration = a;
        if(condition == true){
            current_planet.K = kinetic_energy(current_planet);
            current_planet.moment = single_momentum(current_planet);
            current_planet.L = angular_momentum(current_planet);

        }
    }
    condition = false;
}

double quantities::beta_force(body current_body, body other_body,double distance, double beta){
    double f;
    f = -(g*current_body.mass*other_body.mass)/(pow(distance,beta));
    return f;
}

double quantities::kinetic_energy(body current_body){

    double Kin = 0.5*current_body.mass*current_body.velocity.lengthSquared();
    return Kin;
}

void quantities::COM(){
    R_CoM.zeros();
    double total_mass = 0;
    vec3 r;



}



double quantities::single_momentum(body current_body){
    double v = current_body.velocity.length();
    return current_body.mass*v;
}

void quantities::calculate_potential_energy(){

    for(int i = 0; i < list_of_bodies.size(); i++){
        body &current_planet = list_of_bodies[i];
        total_potential_energy = 0;
        for(int j = i+1; j < list_of_bodies.size(); j++){
            body &other_planet = list_of_bodies[j];
            vec3 relative_coordinates = current_planet.position - other_planet.position;
            total_potential_energy  += -(g*current_planet.mass*other_planet.mass)/(relative_coordinates.length());

        }
        current_planet.P = total_potential_energy;
    }
}

double quantities::angular_momentum(body current_body){
    vec3 l = (current_body.position).cross(current_body.mass*current_body.velocity);
    double ll = l.length();
    return ll;
}
void quantities::write_to_file(double beta){
    for(store_values info: information){
        info.write_to_file(beta);
    }
}

double quantities::relativistic_force(body current_planet, body other_planet, double distance){
    double c = 63239.7263;
    double rel_force = (-g*current_planet.mass*other_planet.mass/(distance*distance))*(
                1 + (3*current_planet.L*current_planet.L)/(distance*distance*c*c));
    return rel_force;
}



