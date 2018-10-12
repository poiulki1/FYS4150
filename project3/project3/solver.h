#ifndef SOLVER_H
#define SOLVER_H
#include <body.h>

#define _USE_MATH_DEFINES
#include <cmath>

class solver
{
public:
    friend class body;

    vector<body> list_of_bodies;
    int total_bodies;
    double g_const;

    solver();

    //functions
    void add_new_body(body new_body);
    void euler_one_step(double dt);

    // hentet fra body klassen
    double distance(body other_body);
    vector<vec3> acceleration(double g_constant);
    double kinEnergy();
    double potEnergy(body &other_body, double g_constant, double eps);
};

#endif // SOLVER_H
