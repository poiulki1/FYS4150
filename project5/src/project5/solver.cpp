#include "solver.h"

solver::solver(int Mc, double initial_S, double initial_I,double initial_R,
               int N_population,double final_time, double start_time):
    MCcycles(Mc),
    N(N_population),
    ft(final_time),
    ts(start_time),

    init_S(initial_S),
    init_I(initial_I),
    init_R(initial_R)
{

    parameter = new double[7]();

}

solver::~solver(){
    delete[] parameter;

}


void solver::initialize_parameters( double a, double b, double c,
                                    double d, double d_i, double e, double f,
                                    double deviation, double w){
    parameter[0] = a;
    parameter[1] = b;
    parameter[2] = c;
    parameter[3] = d;
    parameter[4] = d_i;
    parameter[5] = e;
    parameter[6] = f;
    a0 = a;
    f0 = f;

    limit = deviation;
    omega = w;
}

void solver::MonteCarlo(double dt){

    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    //moving from S to I to R
    S_to_I = (parameter[0]*S*I/N)*dt;
    I_to_R = parameter[1]*I*dt;
    R_to_S = parameter[2]*R*dt;
    S_to_R = parameter[6]*dt; //vaccination

    //death which occours in group S,I,R

    S_to_D = (parameter[3]*S)*dt;
    I_to_Di = (parameter[4])*dt*I;
    I_to_D =  (parameter[3])*I*dt;
    R_to_D = (parameter[3]*R)*dt;

    //babies born into group S
    B_to_S = (parameter[5]*dt*N);



    susceptible = infected = recovered = 0;

    if(RandomNumberGenerator(gen) < S_to_I){
        susceptible -= 1;
        infected += 1;
    }

    if(RandomNumberGenerator(gen) < I_to_R){
        infected -= 1;
        recovered += 1;
    }

    if(RandomNumberGenerator(gen) < R_to_S){
        recovered -= 1;
        susceptible += 1;

    }

    if(RandomNumberGenerator(gen) < S_to_D){
        susceptible -= 1;
    }

    if(RandomNumberGenerator(gen) < I_to_Di){
        infected -= 1;

    }

    if(RandomNumberGenerator(gen) < I_to_D){
        infected -= 1;
    }

    if(RandomNumberGenerator(gen) < R_to_D){
        recovered -= 1;

    }

    if(RandomNumberGenerator(gen) < B_to_S ){
        susceptible += 1;

    }

    if(RandomNumberGenerator(gen) < S_to_R){
        susceptible -= 1;
        recovered += 1;

    }

    S += susceptible;
    I += infected;
    R += recovered;

    Npopulation(S, I, R);

}

rowvec solver::derivative(int time, rowvec &y){
    rowvec pd_f_t = zeros<rowvec>(3);

    s = y(0);
    i = y(1);
    r = y(2);

    pd_f_t(0) = parameter[2]*r - parameter[0]*(s*i)/N - parameter[3]*s + parameter[5]*(s+i+r) - parameter[6];
    pd_f_t(1) = parameter[0]*(s*i)/N - parameter[1]*i - parameter[3]*i - parameter[4]*i;
    pd_f_t(2) = parameter[1]*i - parameter[2]*r - parameter[3]*r + parameter[6];
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

    return dt_step;

}

void solver::update_a(double time){
    parameter[0] = limit*cos(omega*time) + a0;
    //cout <<parameter[0] << "----" << limit*cos(omega*time) <<"----" <<time <<"----"  << N <<endl;
}

void solver::update_f(double time){
    parameter[6] = f0 + f0*log(1+time);

    //cout <<parameter[0] << "----" << limit*cos(omega*time) <<"----" <<time <<"----"  << N <<endl;
}

void solver::change_b(double value){
    parameter[1] = value;
}

void solver::Npopulation(double S_, double I_, double R_){
    N = (S_ + I_ + R_);


}

void solver::set_initial(bool Mc_arg, mat &M){
    if(Mc_arg == true){
        S = init_S;
        I = init_I;
        R = init_R;
    }
    else{

        M(0,0) = init_S;
        M(0,1) = init_I;
        M(0,2) = init_R;

    }
}

void solver::reset_matrix(mat &A, rowvec &B){
    A.zeros();
    B.zeros();
}

void solver::calculate_variance_and_sigma(int b, rowvec &time, mat &sampled_SIR, mat &unsampled_SIR,
                                          int nsamples, savetofile &save_obj,int length){

    mat variance = zeros<mat>(MCcycles+1, number_of_functions);
    rowvec sigma = zeros<rowvec>(number_of_functions);
    string arg = "variance";

    for(int j = 0; j < length; j++){
        variance(j,0) += (sampled_SIR(j,0)-unsampled_SIR(j,0))*(sampled_SIR(j,0)-unsampled_SIR(j,0))/(nsamples-1);
        variance(j,1) += (sampled_SIR(j,1)-unsampled_SIR(j,1))*(sampled_SIR(j,1)-unsampled_SIR(j,1))/(nsamples-1);
        variance(j,2) += (sampled_SIR(j,2)-unsampled_SIR(j,2))*(sampled_SIR(j,2)-unsampled_SIR(j,2))/(nsamples-1);

        sigma(0) += sqrt(variance(j,0))/(nsamples*length);
        sigma(1) += sqrt(variance(j,1))/(nsamples*length);
        sigma(2) += sqrt(variance(j,2))/(nsamples*length);
    }
    save_obj.save(arg, b, variance,time, length );
    cout << "The standard diviation in S = " <<setprecision(5) <<sigma(0) << "|" << "I = "
         <<setprecision(5) << sigma(1) << "|" << "R = "
           <<setprecision(5) <<sigma(2) << endl;
    //save_obj.save_standard_diviation(sigma(0), sigma(1), sigma(2), b);
}

void solver::execute_solve(string arg, bool Mc_arg, double B, int nsamples, savetofile &save_obj){
    double h;
    if(Mc_arg == true){
        MCcycles = MCcycles*10;
        h = double(ft-ts)/MCcycles;
    }
    else{
        MCcycles = MCcycles;
        h = double(ft-ts)/MCcycles;
    }

    rowvec time = zeros<rowvec>(MCcycles+1);
    rowvec step_array = zeros<rowvec>(nsamples);
    mat SIR = zeros<mat>(MCcycles+1, number_of_functions);
    mat raw_SIR = zeros<mat>(MCcycles+1, number_of_functions);



    for(int b = 1; b < B; b++){
        reset_matrix(SIR, time);
        int length = 0;
        if(Mc_arg == true){;
            for(int n = 0; n < nsamples; n++){

                set_initial(Mc_arg, SIR);
                Npopulation(S,I,R);
                change_b(b);
                int mc_step = 0;
                while(time(mc_step) <= ft + 1){
                    update_a(time(mc_step));
                    update_f(time(mc_step));

                    SIR(mc_step,0) += S/double(nsamples);
                    SIR(mc_step,1) += I/double(nsamples);
                    SIR(mc_step,2) += R/double(nsamples);



                    raw_SIR(mc_step,0) =  S;
                    raw_SIR(mc_step,1) =  I;
                    raw_SIR(mc_step,2) =  R;

                    time(mc_step+1) = time(mc_step) + dt(N);

                    MonteCarlo(dt(N));


                    if(N <= 1){
                        MCcycles = mc_step;
                        break;
                    }
                    mc_step += 1;

                }
                step_array(n) = mc_step;

            }

            save_obj.save(arg, b, SIR,time, step_array.min());
            calculate_variance_and_sigma(b, time, SIR, raw_SIR,nsamples, save_obj, step_array.min());
        }

        else{

            change_b(b);
            reset_matrix(SIR, time);
            set_initial(Mc_arg, SIR);
            Npopulation(SIR(0,0),SIR(0,1),SIR(0,2));


            double h = double(ft-ts)/MCcycles;
            for(int i = 0; i < MCcycles; i++){



                update_a(time(i));
                update_f(time(i));

                rk4(time(i), SIR, 0, i, h);
                rk4(time(i), SIR, 1, i, h);
                rk4(time(i), SIR, 2, i, h);

                Npopulation(SIR(i,0),SIR(i,1),SIR(i,2));
                time(i + 1) = time(i) + h;
            }
            save_obj.save(arg, b, SIR,time,MCcycles);
        }
    }
}


