#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <armadillo>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace arma;
using namespace std;

inline int periodic(int i, int limit){return (i+limit) % (limit);}

void initialize(int, double, mat&, double&, double&);
void metropolis(int, mat&, double&, double&, double *);
void write_to_file(ofstream &, int, int, double, vec);
void init_matrix(mat&, int, string);

#endif // METROPOLIS_H
