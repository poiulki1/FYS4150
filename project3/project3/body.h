#ifndef BODY_H
#define BODY_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;

class body
{
public:
    //Properties
    double mass;
    double position[3];
    double velocity[3];
    double potential_energy;
    double kinetic_energy;

    //Initializers
    body();
    body(double M, double x, double y, double z, double vx, double vy, double vz);

    //Functions
    double distance(body other_body);
    double grav_force(body other_body, double g_constant);
    double acceleration(body other_body, double g_constant);
    double kinEnergy();
    double potEnergy(body &other_body, double g_constant, double eps);
};

#endif // BODY_H
