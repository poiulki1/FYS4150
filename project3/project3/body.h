#ifndef BODY_H
#define BODY_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <vec3.h>
using std::vector;

class body
{
public:
    //Properties
    double mass;
    vec3 position;
    vec3 velocity;
    double potential_energy;
    double kinetic_energy;
    double m_constant;

    //Initializers
    body();
    body(double M, vec3 pos, vec3 vel);
    //body(double M, double x, double y, double z, double vx, double vy, double vz);

    //Functions
    double distance(body other_body);
    /*
    vec3 grav_force(body other_body, double g_constant);
    vec3 acceleration(body other_body, double g_constant);
    double kinEnergy();
    double potEnergy(body &other_body, double g_constant, double eps);*/
};

#endif // BODY_H
