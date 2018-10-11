#include "body.h"
#include <armadillo>

body::body(arma::vec pos, arma::vec vel, double mass)
{
    b_mass = mass;
    b_position = pos;
    b_velocity = vel;
}
