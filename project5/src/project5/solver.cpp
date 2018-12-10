#include "solver.h"

solver::solver(int Mc, double initial_S, double initial_I,double initial_R, double a, double b, double c,
               double d, double d_i, double e, double f, int N_population,double final_time, double start_time):
    MCcycles(Mc),
    N(N_population),
    ft(final_time),
    ts(start_time),

    init_S(initial_S),
    init_I(initial_I),
    init_R(initial_R)
{

    parameter = new double[7]();

    //rowvec kk2 = zeros(number_of_functions);
    //rowvec kk3 = zeros(number_of_functions);
    //rowvec kk4 = zeros(number_of_functions);


    parameter[0] = a;
    parameter[1] = b;
    parameter[2] = c;
    parameter[3] = d;
    parameter[4] = d_i;
    parameter[5] = e;
    parameter[6] = f;


}

solver::~solver(){
    delete[] parameter;

}


void solver::MonteCarlo(int index, double dt, int N, int nsamples){
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    //moving from S to I to R
    S_to_I = (parameter[0]*S*I/N)*dt;
    I_to_R = parameter[1]*I*dt;
    R_to_S = parameter[2]*R*dt;

    //death which occours in group S,I,R
    S_to_D = (parameter[3]*S)*dt;
    I_to_D = (parameter[4] + parameter[3])*I*dt;
    R_to_D = (parameter[3]*R)*dt;

    //babies born into group S
    B_to_S = (dt*parameter[5]*N);


    ss = ii = rr = dd = 0;

    if(RandomNumberGenerator(gen) < S_to_I){
        ss -= 1;
        ii += 1;
    }

    if(RandomNumberGenerator(gen) < I_to_R ){
        ii -= 1;
        rr += 1;
    }

    if(RandomNumberGenerator(gen) < R_to_S){
        rr -= 1;
        ss += 1;

    }

    if(RandomNumberGenerator(gen) < S_to_D){
        ss -= 1;
    }

    if(RandomNumberGenerator(gen) < I_to_D){
        ii -= 1;

    }

    if(RandomNumberGenerator(gen) < R_to_D){
        rr -= 1;

    }

    if(RandomNumberGenerator(gen) < B_to_S){
        ss += 1;

    }

    S += ss;
    I += ii;
    R += rr;

    Npopulation(S, I, R);
    cout<< N << endl;



}

rowvec solver::derivative(int time, rowvec &y){
    rowvec pd_f_t = zeros<rowvec>(3);
    s = y(0);
    i = y(1);
    r = y(2);

    pd_f_t(0) = parameter[2]*r - parameter[0]*(s*i)/N - parameter[3]*s + parameter[5]*i - parameter[6];
    pd_f_t(1) = parameter[0]*(s*i)/N - parameter[1]*i - parameter[3]*i - parameter[4]*i;
    pd_f_t(2) = parameter[1]*i - parameter[2]*r - parameter[3]*i + parameter[6];
    return pd_f_t;
}

void solver::rk4(int time, mat &y, int dim, int index,double h){

    rowvec K0_fac = zeros<rowvec>(number_of_functions);
    rowvec K1_fac = zeros<rowvec>(number_of_functions);
    rowvec K2_fac = zeros<rowvec>(number_of_functions);
    rowvec K3_fac = zeros<rowvec>(number_of_functions);
    mat K = zeros<mat>(4,number_of_functions);


    K0_fac = y.row(index);
    K.row(0) = h*derivative(time, K0_fac);

    K1_fac = K.row(0)*0.5 + y.row(index);

    K.row(1) = h*derivative(time + h/2, K1_fac);

    K2_fac = K.row(1)*0.5 + y.row(index);
    K.row(2) = h*derivative(time + h/2, K2_fac);

    K3_fac = K.row(2) + y.row(index);
    K.row(3) = h*derivative(time + h, K3_fac);

    y(index + 1, dim) = y(index,dim) + (K(0,dim) + 2*K(1,dim) + 2*K(2,dim) + K(3,dim))/6;

}

double solver::dt(double N){
    double dt_step;
    rowvec dt_array = zeros<rowvec>(7);
    dt_array(0) = dt1(parameter[0],N);
    dt_array(1) = dt2(parameter[1],N);
    dt_array(2) = dt3(parameter[2],N);
    dt_array(3) = dt4(parameter[3],N);
    dt_array(4) = dt5(parameter[4],parameter[3] ,N);
    dt_array(5) = dt6(parameter[3],N);
    dt_array(6) = dt7(parameter[5],N);

    dt_step = dt_array.min();


    /*
    double dt_step;
    if(dt2(value2, N) <= dt1(value1, N) && dt2(value1, N) < dt3(value3, N)){
        dt_step =  dt2(value2, N);
    }
    else if( dt1(value1, N) < dt2(value2, value4, value5, N) && dt1(value1, N) < dt3(value3, N)){
        dt_step = dt1(value1, value4, N);
    }
    else{
        dt_step = dt3(value3, N);
    }
    */

    return dt_step;

}
void solver::change_b(double value){
    parameter[1] = value;
}

void solver::Npopulation(double S_, double I_, double R_){
    N = (S_ + I_ + R_);
}
void solver::execute_solve(bool Mc_arg, double B){

    rowvec tt = zeros<rowvec>(MCcycles+1);
    mat SIR = zeros<mat>(MCcycles+1, number_of_functions);

    int nsamples = 10;

    for(int b = 1; b < B; b++){
        ofstream file;
        file.open("test_" + to_string(b) + ".txt");
        if(Mc_arg == true){
            for(int n = 0; n < nsamples; n++){
                S = init_S;
                I = init_I;
                R = init_R;
                Npopulation(S,I,R);

                for(int mc_step = 0; mc_step < MCcycles; mc_step++){
                    SIR(mc_step,0) += S/double(nsamples);
                    SIR(mc_step,1) += I/double(nsamples);
                    SIR(mc_step,2) += R/double(nsamples);
                    tt(mc_step+1) = tt(mc_step) + dt(N);


                    MonteCarlo(mc_step, dt(N),  N, nsamples);

                }
            }
            for(int save_step = 0; save_step < MCcycles; save_step++){
                file <<setw(12) <<tt(save_step) << setw(12) <<  SIR(save_step, 0) <<
                       setw(12) <<  SIR(save_step, 1) <<
                       setw(12) <<  SIR(save_step, 2) << endl;
            }
            file.close();
        }

        else{
            change_b(b);
            SIR(0,0) = init_S;
            SIR(0,1) = init_I;
            SIR(0,2) = init_R;

            double h = double(1.0/10000.0);
            for(int i = 0; i < MCcycles; i++){
                rk4(i, SIR, 0, i, h);
                rk4(i, SIR, 1, i, h);
                rk4(i, SIR, 2, i, h);
                tt(i + 1) = tt(i) + h;


                file << setw(12) << i << setw(12)<< SIR(i,0)
                     << setw(12) << SIR(i,1) << setw(12) << SIR(i,2) << endl;

            }

        }


    }
}


/*
 * void solver::derivatives(double t, double s, double i , double r){

    sir[0] = parameter[2]*r - parameter[0]*(s*i)/N - parameter[3]*s + parameter[5]*i - parameter[6];

    sir[1] = parameter[0]*(s*i)/N - parameter[1]*i - parameter[3]*i - parameter[4]*i;

    sir[2] = parameter[1]*i - parameter[2]*r - parameter[3]*i + parameter[6];


}


void solver::rk4_step(double h, double time,int dim, int index){
    double sir[0] = y[idx(0, index)];
    double sir[1] = y[idx(1, index)];
    double sir[2] = y[idx(2, index)];

    double s = y[idx(0, index)];
    double i = y[idx(1, index)];
    double r = y[idx(2, index)];

    for(int i = 0; i < 3; i++){
        k1 = derivatives(time, s, i, r); //k1 - updates the sir vector

        k2 = derivatives(time + h*0.5, s + h*sir[0]*0.5, i + h*sir[1]*0.5, r  + h*sir[2]*0.5);

        k3 = derivatives(time + h*0.5, s  + h*sir[0]*0.5, i  + h*sir[1]*0.5, r + h*sir[2]*0.5);

        k4 = derivatives(time + h  , s  + h*sir[0], i + h*sir[1],  + h*sir[2]);
    }

    y[idx(dim, index + 1)] = y[idx(dim, index)] + (h/6)*(k1 + 2*k2 + 2*k3 + k4); //1/6 eller h/6 ????? hvis h'ene er på
                                                                                 //i k uttrykene trenger vi ikke det her
}
        else{
            change_b(b);
            cout << " Initial S " << y[idx(0,0)] << " Initial I " << y[idx(1,0)]
                 << " Initial R " << y[idx(2,0)] << endl;
            //cout << parameter[1] << " " << b << endl;
            for(int mc_step = 0; mc_step < MCcycles; mc_step++){



                rk4_step(0.00001, mc_step, 0, mc_step);
                rk4_step(0.00001, mc_step, 1, mc_step);
                rk4_step(0.00001, mc_step, 2, mc_step);

                file << setw(12) << mc_step << setw(12)<< y[idx(0,mc_step)]
                     << setw(12) << y[idx(1,mc_step)] << setw(12) << y[idx(2,mc_step)] << endl;
            }
            file.close();
            cout << " end S " << y[idx(0,MCcycles)] << " end I " << y[idx(1,MCcycles)]
                 << " end R " << y[idx(2,MCcycles)]
                 << " total " << y[idx(0,MCcycles)] + y[idx(1,MCcycles)] + y[idx(2,MCcycles)] << endl;

            for(int restet_array = 0; restet_array < MCcycles; restet_array++){
                y[idx(0,restet_array)] = 0;
                y[idx(1,restet_array)] = 0;
                y[idx(2,restet_array)] = 0;
            }



            y[idx(0,0)] = 300;
            y[idx(1,0)] = 100;
            y[idx(2,0)] = 0;
        }

    }

void solver::change_b(double value){
    parameter[1] = value;
}

void solver::execute_solve(bool Mc_arg, double B){
    double a = parameter[0];
    double c = parameter[2];


    mat SIR = zeros<mat>(MCcycles+1, number_of_functions);
    SIR(0,0) = 300;
    SIR(0,1) = 100;
    SIR(0,2) = 0;
    double time = double((ft - ts))/MCcycles;
    for(int b = 1; b < B; b++){
        ofstream file;
        file.open("test_" + to_string(b) + ".txt");

        if(Mc_arg == true){
            for(int mc_step = 0; mc_step < MCcycles; mc_step++){

                if(dt3(b) <= dt1(a) && dt3(b) < dt2(c)){

                    S_to_I[mc_step] = (a*S[mc_step]*I[mc_step]/N)*dt3(b);
                    I_to_R[mc_step] = b*I[mc_step]*dt3(b);
                    R_to_S[mc_step] = c*R[mc_step]*dt3(b);
                    tid[mc_step+1] = tid[mc_step] + dt3(b);



                }
                else if( dt1(a) < dt3(b) && dt1(a) < dt2(c)){

                    S_to_I[mc_step] = (a*S[mc_step]*I[mc_step]/N)*dt1(a);
                    I_to_R[mc_step] = b*I[mc_step]*dt1(a);
                    R_to_S[mc_step] = c*R[mc_step]*dt1(a);
                    tid[mc_step+1] = tid[mc_step] + dt1(a);




                }
                else{

                    S_to_I[mc_step] = (a*S[mc_step]*I[mc_step]/N)*dt2(c);
                    I_to_R[mc_step] = b*I[mc_step]*dt2(c);
                    R_to_S[mc_step] = c*R[mc_step]*dt2(c);
                    tid[mc_step+1] = tid[mc_step] + dt2(c);

                }
                cout << S_to_I[mc_step] << " " << I_to_R[mc_step] << " " << R_to_S[mc_step] << endl;
                MonteCarlo(mc_step,  S_to_I[mc_step],  I_to_R[mc_step] ,  R_to_S[mc_step], N);
                file << setw(12) << tid[mc_step] << setw(12)<< S[mc_step]
                        << setw(12) << I[mc_step] << setw(12) << R[mc_step] << endl;
            }
            file.close();
        }
        else{
            change_b(b);
            for(int i = 0; i < MCcycles; i++){
                rk4(i, SIR, 0, i);
                rk4(i, SIR, 1, i);
                rk4(i, SIR, 2, i);

                file << setw(12) << i << setw(12)<< SIR(i,0)
                     << setw(12) << SIR(i,1) << setw(12) << SIR(i,2) << endl;

            }

        }


    }
                    S_to_I[mc_step] = (a*S[mc_step]*I[mc_step]/N)*dt2(b);
                    I_to_R[mc_step] = b*I[mc_step]*dt2(b);
                    R_to_S[mc_step] = c*R[mc_step]*dt2(b);
                    tid[mc_step+1] = tid[mc_step] + dt2(b);
                }
                else if( dt1(a) < dt2(b) && dt1(a) < dt3(c)){
                    S_to_I[mc_step] = (a*S[mc_step]*I[mc_step]/N)*dt1(a);
                    I_to_R[mc_step] = b*I[mc_step]*dt1(a);
                    R_to_S[mc_step] = c*R[mc_step]*dt1(a);
                    tid[mc_step+1] = tid[mc_step] + dt1(a);
                }
                else{
                    S_to_I[mc_step] = (a*S[mc_step]*I[mc_step]/N)*dt3(c);
                    I_to_R[mc_step] = b*I[mc_step]*dt3(c);
                    R_to_S[mc_step] = c*R[mc_step]*dt3(c);
                    tid[mc_step+1] = tid[mc_step] + dt3(c);
                }
}
*/



