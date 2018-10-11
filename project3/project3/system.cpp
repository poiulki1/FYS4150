#include "system.h"
#include <armadillo>
#include "body.h"

system::system()
{
    void system::init_bodies(arma::vec initial_pos, arma::vec initial_vel, double mass){
        list_of_bodies.push_back(body(initial_pos, initial_vel, mass));
        num_of_bodies = list_of_bodies.size();
    }
}
