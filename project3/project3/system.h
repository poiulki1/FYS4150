#ifndef SYSTEM_H
#define SYSTEM_H
#include <armadillo>
#include <math.h>
#include <vector>



class system
{
public:
    system();

    void init_bodies();
    int num_of_bodies;
    void find_acceleration();
    double G = 4*M_PI*M_PI;

    list_of_bodies<body> = list_of_bodies;

};

#endif // SYSTEM_H
