#include <iostream>
using namespace std;
#include "solver.h"
#include "vec3.h"
#include "body.h"

solver::solver()
{
    g_const = (4*M_PI*M_PI);
}

void solver::add_new_body(body new_body)
{
    list_of_bodies.push_back(new_body);
    //cout << list_of_bodies.size() << endl;
}

void solver::euler_one_step(double dt)
{
    for(int body_nr=0; body_nr<list_of_bodies.size(); body_nr++){

        //cout << "Planet nr. " << body_nr << endl;


        body &current = list_of_bodies[body_nr];

        vector<vec3> acc = acceleration(g_const);

        //cout << "Acceleration" << acc[body_nr] << endl;

        //integration, one step of Euler Forward
        current.position += current.velocity*dt;
        current.velocity += acc[body_nr]*dt;

        //cout << "Position: " << current.position << endl;
        //cout << "Velocity: " << current.velocity << endl;

        list_of_bodies[body_nr] = current;
    }
}

vector<vec3> solver::acceleration(double g_constant){
    // calculate gravitational force between two bodies as a vector
    int num_bodies = list_of_bodies.size();

    vector<vec3> acc(num_bodies);
    vector<vec3> force(num_bodies);

    for(int i_body = 0; i_body < num_bodies; i_body++){
        cout << "first loop planet nr: " << i_body << endl;
        body body1 = list_of_bodies[i_body];
        for(int j_body = i_body+1; j_body < num_bodies; j_body++){
            cout << "second loop planet nr: " << j_body << endl;
            body body2 = list_of_bodies[j_body];

            vec3 rel_coordinates = body1.position - body2.position;
            double length_of_rel_coord = rel_coordinates.length();
            vec3 rel_coord_unit = rel_coordinates/length_of_rel_coord;

            double rel_distance = body1.distance(body2);

            if(rel_distance != 0){
                double weight = (g_constant*body1.mass*body2.mass)/(rel_distance*rel_distance);
                force[i_body] += (weight*rel_coord_unit);
                acc[i_body] += (weight*rel_coord_unit)/body1.mass;
            }
            else{
                cout << "NULL acc for planet nr.: " << i_body << endl;
                acc[i_body] = vec3(0.,0.,0.);
                force[i_body] = vec3(0.,0.,0.);
            }
        }
    }
    return acc;
}

/*
double solver::kinEnergy(){
    double vel2 = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]);
    return 0.5*this->mass*vel2;
}

double solver::potEnergy(body &other_body, double g_constant, double eps){
    if(eps == 0.0){
        return -(g_constant*this->mass*other_body.mass)/this->distance(other_body);
    }
    else{
        return (g_constant*this->mass*other_body.mass/eps)*(atan(this->distance(other_body)/eps) - (0.5*M_PI));
    }
}*/
