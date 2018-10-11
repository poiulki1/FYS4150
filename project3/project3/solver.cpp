#include "solver.h"
#include "vec3.h"

solver::solver()
{
    void Euler_solver(int N, double dt, vec3 pos,vec3 velocity, vec3 time, vec3 acceleration){
        for(int i = 0; i<=N; i++){
            position(i+1) = position(i) + velocity(i)*dt;

            velocity(i+1) = velocity(i) + acceleration(position(i))*dt;
            time(i+1) = time(i) + dt;
        }
    }

    void Velocity_verlet(int N, double dt,vec3 position, vec3 velocity, vec3 time, vec3 acceleration){
        for(int j = 0; j <= N; j++){
            position(j+1) = position(j) + velocity(j)*dt + 0.5*acceleration(j)*dt*dt;
            velocity(j+1) = velocity(j) + 0.5*acceleration(i)*dt + 0.5*acceleration(j+1)*dt;
            time(j+1) = time(j) + dt;
        }
    }
}
