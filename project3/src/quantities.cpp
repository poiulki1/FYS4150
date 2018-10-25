#include "quantities.h"
#include <iostream>
#include <iomanip>

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

        if(condition == true){
            current_planet.K = kinetic_energy(current_planet);
            current_planet.moment = single_momentum(current_planet);
            current_planet.L = angular_momentum(current_planet);
        }




        for(int j = 0; j < list_of_bodies.size(); j++){
            if(i == j){
                continue;
            }

            body &other_planet = list_of_bodies[j];
            vec3 relative_coordinates = current_planet.position - other_planet.position;
            vec3 unit_vector = relative_coordinates/relative_coordinates.length();
            distance = relative_coordinates.length();
            //current_planet.angle = calculate_perihelion_angle(current_planet, relative_coordinates);


            a += force(current_planet, other_planet, distance)*unit_vector/current_planet.mass;
            //a += beta_force(current_planet, other_planet,distance, beta)*unit_vector/current_planet.mass;

            //a += relativistic_force(current_planet, other_planet, distance)*unit_vector/current_planet.mass;




        }

        current_planet.acceleration = a;




    }
    condition = true;
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



double quantities::calculate_perihelion_angle(body current_body, vec3 R){

    prev_prev_pos = prev_pos;
    prev_prev_vec = prev_vec;

    prev_pos = current_pos;
    prev_vec = vecR;

    current_pos = R.length();
    vecR = R;

    cout << prev_prev_pos << " | " << prev_pos << " | " << current_pos << endl;

    if(prev_pos < prev_prev_pos && prev_pos < current_pos){

        return atan2(prev_vec[1],prev_vec[0]);

    }
    else{
        return 0;
    }
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

    double r_x_v = ((current_planet.position).cross(current_planet.velocity)).lengthSquared();

    double newtonian_force = (-g*current_planet.mass*other_planet.mass/(distance*distance));
    double peturbation_const = (3*r_x_v)/(distance*distance*c*c);

    double rel_force = peturbation_const*newtonian_force;

    //cout << std::setprecision(12) <<peturbation_const <<endl;//<< " " << newtonian_force << endl;
    return rel_force;
}



