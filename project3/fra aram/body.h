#ifndef BODY_H
#define BODY_H
#include <cmath>
#include <vector>
#include <vec3.h>
using namespace std;

class body
{
public:

    vec3 position;
    vec3 velocity;
    vec3 acceleration;
    vec3 prev_acceleration;
    double P;
    double K;
    double L;
    double moment;
    double angle;
    double mass;
    string argument;


    body(double m, vec3 pos, vec3 vel, string arg);

};

#endif // BODY_H
