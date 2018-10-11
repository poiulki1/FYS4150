#ifndef BODY_H
#define BODY_H


class body
{
public:
    body(arma::vec pos, arma::vec vel, double mass);
    double b_mass;
    arma::vec b_position;
    arma::vec b_velocity;
    arma::vec acceleration = arma::vec<arma::zeros>(2);
};

#endif // BODY_H
