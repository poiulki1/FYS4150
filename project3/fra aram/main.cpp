#include <iostream>
#include "vec3.h"
#include "body.h"
#include "solver.h"
#include "quantities.h"
#include "time.h"
#include <fstream>
using namespace std;

int main()
{
    /*
    int esc_vel = 3;
    bool rel_force = true;

    if(esc_vel == 1){
        vector<double>v = {2*M_PI,7.0, 7.50,7.90,8.20,8.50,8.89};
        for(int k = 0; k < v.size(); k++)
        {
            vec3 earth_pos(1.0,0.0,0.0);
            vec3 earth_vel(0.0,v[k],0.0);

            vec3 sun_pos(0.0,0.0,0.0);
            vec3 sun_vel(0.0,0.0,0.0);

            body Earth(0.000003, earth_pos, earth_vel, "earth");
            body Sun(1.0, sun_pos, sun_vel, "sun");

            int N = 10000;
            double final_time = 14.0;
            //string arg = "foward_euler";
            string arg = "velocity_verlet";
            double dt = final_time/((double) N);

            quantities solar_system;
            solar_system.add_new_body(Sun, store_values(N));
            solar_system.add_new_body(Earth, store_values(N));
            solver solutions(dt, solar_system);

            double beta = 2.0;
            string condition2 = "reg_force";
            for(int i=0; i<N; i++){
                solutions.velocity_verlet(solar_system, i, beta);
            }
            solar_system.write_to_file(v[k]);
        }
    }
    if(esc_vel == 2){
        vec3 saturn_pos(1.569973482245909,-9.931807043056629,1.101934926120388e-1);
        vec3 saturn_vel(5.202769599394199e-3,8.527763066533207e-4,-2.216361381984069e-4);
        saturn_vel = saturn_vel*365.25;

        vec3 mars_pos(1.381399662147958E+00,-1.174156225640810E-01,-3.658517548407784E-02);
        vec3 mars_vel(1.784568140790161E-03,1.513529311587753E-02, 2.733079336331583E-04);

        mars_vel = mars_vel*365.25;



        vec3 jupiter_pos(-2.647556619196740E+00, -4.665803344838063E+00,7.857338643883652e-02);
        vec3 jupiter_vel(6.472894738201937e-03, -3.365002410907006e-03, -1.308549331239642e-04);

        vec3 earth_pos(9.221407681304493e-01, 3.859513377832802e-01, -9.342283541098679e-05);
        vec3 earth_vel(-6.828707544910924e-03, 1.585182505936633e-02, 8.139769320198632e-08);



        earth_vel = earth_vel*365.25;
        jupiter_vel = jupiter_vel*365.25;



        vec3 sun_pos(-1.626563776599144E-04, 7.254085928335599E-03, -7.230182567169775E-05);
        vec3 sun_vel(-7.590756043911132E-06,2.572902340170582E-06,1.897452701696540E-07);
        sun_vel = sun_vel*365.25;

        body Earth(0.000003, earth_pos, earth_vel, "earth");
        body Sun(1.0, sun_pos, sun_vel, "sun");
        body Jupiter(0.000949, jupiter_pos,jupiter_vel, "jupiter");
        body Mars(3.27e-7, mars_pos, mars_vel, "mars");
        body Saturn(0.00028517, saturn_pos, saturn_vel, "saturn");

        int N = 10000;
        double final_time = 4;
        //string arg = "foward_euler";
        string arg = "velocity_verlet";
        double dt = final_time/((double) N);

        quantities solar_system;
        solar_system.add_new_body(Sun, store_values(N));
        solar_system.add_new_body(Earth, store_values(N));

        //solar_system.add_new_body(Jupiter,store_values(N));
        //solar_system.add_new_body(Mars,store_values(N));
        //solar_system.add_new_body(Saturn, store_values(N));



        solver solutions(dt, solar_system); //Ta inn solar system i konstruktoren



        cout << "earth :" <<solar_system.list_of_bodies[1].position << endl;
        //cout << "earth :" <<solar_system.list_of_bodies[1].position << endl;
        //cout << "jupiter :"<<solar_system.list_of_bodies[2].position << endl;
        //cout << "saturn :"<<solar_system.list_of_bodies[3].position << endl;


        clock_t start_time, end_time;
        double time_used;

        //string condition2 = "reg_force";
        string condition2 = "beta_force";



        start_time = clock();
        if(condition2 == "reg_force"){
            double beta = 2;
            for(int i=0; i<N; i++){

                solutions.velocity_verlet(solar_system, i, beta);

            }
            end_time = clock();
            solar_system.write_to_file(beta);

            time_used = ((double)(end_time - start_time))/CLOCKS_PER_SEC;
            cout << "time used for N = " << N << "and dt = " << dt << "is t = " << time_used << endl;
        }


        if(condition2 =="beta_force"){
            vector<double> beta = {2.0,2.3,2.5,2.8,3.0};
            for(int j = 0; j < beta.size(); j ++){
                for(int i=0; i<N; i++){

                    solutions.velocity_verlet(solar_system, i, beta[j]);

                }
                solar_system.write_to_file(beta[j]);
            }
        }
        cout << "---------------------------------------------------"<<endl;
        cout << "earth :" <<solar_system.list_of_bodies[1].position << endl;
        //cout << "earth :" <<solar_system.list_of_bodies[1].position << endl;
        //cout << "jupiter :"<<solar_system.list_of_bodies[2].position << endl;
        //cout << "saturn :"<<solar_system.list_of_bodies[3].position << endl;
    }

    if(rel_force == true){
        */


    //vec3 mercury_pos(-1.108553727766392E-01, -4.453177109185471E-01, -2.689875721611897E-02);
    //vec3 mercury_vel(2.167391218861134E-02, -5.267049408145004E-03, -2.419432184542312E-03);
    //mercury_vel = mercury_vel*365.25;

    //vec3 sun_pos(-1.626563776599144E-04, 7.254085928335599E-03, -7.230182567169775E-05);
    //vec3 sun_vel(-7.590756043911132E-06,2.572902340170582E-06,1.897452701696540E-07);

    vec3 jupiter_pos(-2.647556619196740E+00, -4.665803344838063E+00,7.857338643883652e-02);
    vec3 jupiter_vel(6.472894738201937e-03, -3.365002410907006e-03, -1.308549331239642e-04);

    vec3 earth_pos(9.221407681304493e-01, 3.859513377832802e-01, -9.342283541098679e-05);
    vec3 earth_vel(-6.828707544910924e-03, 1.585182505936633e-02, 8.139769320198632e-08);

    vec3 sun_pos(0,0,0);
    vec3 sun_vel(0,0,0);

    earth_vel = earth_vel*365.25;
    jupiter_vel = jupiter_vel*365.25;

    body Earth(0.000003, earth_pos, earth_vel, "earth");
    body Sun(1.0, sun_pos, sun_vel, "sun");
    double Jupiter_mass = 0.000949;
    vector<double> jupiter_mass_list = {Jupiter_mass, Jupiter_mass*10, Jupiter_mass*1000};

    for(int k = 0; k < jupiter_mass_list.size(); k++){


        body Jupiter(jupiter_mass_list[k], jupiter_pos,jupiter_vel, "jupiter");


        //body Mercury(0.0000001651, mercury_pos, mercury_vel, "mercury");

        int N = 1e5;
        double final_time = 100;


        double dt = final_time/((double) N);

        quantities solar_system;
        solar_system.add_new_body(Sun, store_values(N));
        solar_system.add_new_body(Jupiter, store_values(N));
        solar_system.add_new_body(Earth, store_values(N));

        solver solutions(dt, solar_system);

        double beta = 2.0;


        for(int i=0; i<N; i++){
            solutions.velocity_verlet(solar_system, i, beta);

        }

        solar_system.write_to_file(beta);
    }

    return 0;
}
