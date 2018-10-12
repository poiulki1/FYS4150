#include "body.h"
#include "solver.h"
#include <iostream>
using namespace std;

body::body()
{

}

body::body(double M, vec3 pos, vec3 vel)
{
    mass = M;

    position = pos;
    velocity = vel;

    potential_energy = 0.;
    kinetic_energy = 0.;
}

double body::distance(body other_body)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = this->position[0];
    y1 = this->position[1];
    z1 = this->position[2];

    x2 = other_body.position[0];
    y2 = other_body.position[1];
    z2 = other_body.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
}
/*
vec3 body::grav_force(body other_body, double g_constant){
    // calculate gravitational force between two bodies as a vector
    int K = list_of_bodies.size();// (BUG)

    for(int i_body = 0; i_body < K; i_body++){
        body1 = ;//list_of_bodies[i_body];
        for(int j_body = i_body+1; j_body < K; j_body++){
            body2 = ;//list_of_bodies[j_body];

            vec3 rel_coordinates = body1.position - body2.position;
            double rel_distance = body1.distance(body2);

        }
    }

    double r = this->distance(other_body);



    /*
    // calculates grav force for an object with fixed sun in the origin
    if(r != 0.0){
        double weight = ((g_constant*this->mass*other_body.mass)/(r*r));
        return ((this->position.flip())/(this->position.length()))*weight;
    }
    else{
        return vec3(1.0,0.0,0.0);
    }
}

vec3 body::acceleration(body other_body, double g_constant){
    return (this->grav_force(other_body, g_constant)/(this->mass));
}

double body::kinEnergy(){
    double vel2 = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]);
    return 0.5*this->mass*vel2;
}

double body::potEnergy(body &other_body, double g_constant, double eps){
    if(eps == 0.0){
        return -(g_constant*this->mass*other_body.mass)/this->distance(other_body);
    }
    else{
        return (g_constant*this->mass*other_body.mass/eps)*(atan(this->distance(other_body)/eps) - (0.5*M_PI));
    }
}*/
