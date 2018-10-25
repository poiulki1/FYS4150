#include "body.h"
#include "vec3.h"

body::body(double m, vec3 pos, vec3 vel, string arg)
{
    mass = m;
    position = pos;
    velocity = vel;

    acceleration = {0,0,0};
    prev_acceleration = {0,0,0};
    argument = arg;
    K = 0;
    P = 0;
    L = 0;
    moment = 0;
    angle = 0;


}


