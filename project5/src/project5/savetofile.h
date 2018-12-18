#ifndef SAVETOFILE_H
#define SAVETOFILE_H

#include <armadillo>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace arma;

class savetofile
{
private:
    string filename;
    string number_arg;
    ofstream outfile;

public:
    savetofile(string fname, string number_arg);
    void save(string arg, int b, mat &M, rowvec &y, int run_length);
    void save_each_step(double M1, double M2, double M3, double y);
    void save_standard_diviation(int sigma1, int sigma2, int sigma3, int b);
};

#endif // SAVETOFILE_H
