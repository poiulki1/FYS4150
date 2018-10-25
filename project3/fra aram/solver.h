#ifndef SOLVER_H
#define SOLVER_H
#include "body.h"
#include "vec3.h"
#include "quantities.h"
using namespace std;

class solver
{
public:
    //Solar_system* solarSystem;
    quantities* solarSystem;
    friend class body;
    vector<body> list_of_bodies;
    int total_bodies;
    solver(double ddt, quantities& sSystem);
    vec3 dist;

    double dt;



    void forward_euler(class quantities &system, int iter, double beta);
    void velocity_verlet(quantities &system, int iter, double beta);
};

#endif // SOLVER_H
