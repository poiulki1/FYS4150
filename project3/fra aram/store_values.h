#ifndef STORE_VALUES_H
#define STORE_VALUES_H
#include <vector>
#include "vec3.h"
#include "body.h"


using namespace std;

class store_values
{
public:

    store_values(int N);

    //functions
    void write_to_file(double beta);
    void update(body &noe, int i);

    //vectors
    vector<vec3> position;
    vector<vec3> velocity;
    vector<double> kinetic_energy;
    vector<double> potential_energy;
    vector<double> angular_momentum;
    vector<double> momentum;
    vector<double> theta;


    //constants
    string arg;

};

#endif // STORE_VALUES_H
