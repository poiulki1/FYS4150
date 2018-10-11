#include "body.h"

body::body()
{

}

body::body(double M, double x, double y, double z, double vx, double vy, double vz)
{
    mass = M;

    position[0] = x;
    position[1] = y;
    position[2] = z;

    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;

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

double body::grav_force(body other_body, double g_constant){
    double r = this->distance(other_body);
    if(r != 0.0){
        return (g_constant*this->mass*other_body.mass)/(r*r);
    }
    else{
        return 0;
    }
}

double body::acceleration(body other_body, double g_constant){
    double r = this->distance(other_body);
    if(r != 0.0){
        return this->grav_force(other_body, g_constant)/(this->mass*r);
    }
    else{
        return 0;
    }
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
}
