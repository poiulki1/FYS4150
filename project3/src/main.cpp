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

    int esc_vel = 3;
    bool rel_force = true;
    bool entire_solar = true;

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
    string condition2 = "beta_force";

    if(esc_vel == 2){
        if(condition2 =="beta_force"){
            vector<double> beta = {2.0, 2.3, 2.5, 2.8, 3.0};
            for(int j = 0; j < beta.size(); j ++){
                vec3 saturn_pos(1.569973482245909,-9.931807043056629,1.101934926120388e-1);
                vec3 saturn_vel(5.202769599394199e-3,8.527763066533207e-4,-2.216361381984069e-4);
                saturn_vel = saturn_vel*365.25;

                vec3 mars_pos(1.381399662147958E+00,-1.174156225640810E-01,-3.658517548407784E-02);
                vec3 mars_vel(1.784568140790161E-03,1.513529311587753E-02, 2.733079336331583E-04);

                mars_vel = mars_vel*365.25;



                vec3 jupiter_pos(-2.647556619196740E+00, -4.665803344838063E+00,7.857338643883652e-02);
                vec3 jupiter_vel(6.472894738201937e-03, -3.365002410907006e-03, -1.308549331239642e-04);

                //vec3 earth_pos(9.221407681304493e-01, 3.859513377832802e-01, -9.342283541098679e-05);
                //vec3 earth_vel(-6.828707544910924e-03, 1.585182505936633e-02, 8.139769320198632e-08);



                //earth_vel = earth_vel*365.25;
                //jupiter_vel = jupiter_vel*365.25;



                //vec3 sun_pos(-1.626563776599144E-04, 7.254085928335599E-03, -7.230182567169775E-05);
                //vec3 sun_vel(-7.590756043911132E-06,2.572902340170582E-06,1.897452701696540E-07);
                //sun_vel = sun_vel*365.25;

                vec3 earth_pos(0.98, 0, 0);
                vec3 earth_vel(0, 2*M_PI, 0);

                vec3 sun_pos(0, 0, 0);
                vec3 sun_vel(0,0,0);

                body Earth(0.000003, earth_pos, earth_vel, "earth");
                body Sun(1.0, sun_pos, sun_vel, "sun");
                body Jupiter(0.000949, jupiter_pos,jupiter_vel, "jupiter");
                body Mars(3.27e-7, mars_pos, mars_vel, "mars");
                body Saturn(0.00028517, saturn_pos, saturn_vel, "saturn");

                int N = 10000;
                double final_time = 1;
                //string arg = "foward_euler";
                //string arg = "velocity_verlet";
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


                //clock_t start_time, end_time;
                //double time_used;

                //string condition2 = "reg_force";


                /*
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
        */



                for(int i=0; i<N; i++){

                    solutions.velocity_verlet(solar_system, i, beta[j]);

                }
                solar_system.write_to_file(beta[j]);
            }
        }
        cout << "---------------------------------------------------"<<endl;
        //cout << "earth :" <<solar_system.list_of_bodies[1].position << endl;
        //cout << "earth :" <<solar_system.list_of_bodies[1].position << endl;
        //cout << "jupiter :"<<solar_system.list_of_bodies[2].position << endl;
        //cout << "saturn :"<<solar_system.list_of_bodies[3].position << endl;
    }






    if(entire_solar == true){
        //initial values
        vec3 pluto_pos(1.167424090022203E+01,-3.157815544621814E+01,1.530520743281332E-03);
        vec3 uranus_pos(1.715603176548234E+01,1.002252941244031E+01, -1.849341831335910E-01);
        vec3 venus_pos(6.347210141529195E-01, 3.476452881358778E-01,-3.185781781727027E-02);
        vec3 neptun_pos(2.892839739027679E+01, -7.699294884282640E+00,-5.082059114328670E-01);
        vec3 saturn_pos(1.569973482245909,-9.931807043056629,1.101934926120388e-1);
        vec3 mars_pos(1.391002003207167E+00,-3.309806589020768E-03,-3.420150585256754E-02);
        vec3 jupiter_pos(-2.647556619196740E+00, -4.665803344838063E+00,7.857338643883652e-02);
        //vec3 earth_pos(9.221407681304493e-01, 3.859513377832802e-01, -9.342283541098679e-05);
        //vec3 mercury_pos(6.757570266475262E-02, -4.515846783353264E-01, -4.309982188170661E-02);
        vec3 mercury_pos(0.3075,0,0);
        vec3 earth_pos(1.0,0,0);



        vec3 pluto_vel(3.020174112710275E-03, 4.314350730919463E-04,-9.052934256085976E-04);
        vec3 uranus_vel(-2.006319067923562E-03, 3.209583619757120E-03, 3.784939160901838E-05);
        vec3 venus_vel(-9.780339136705398E-03, 1.765226772287646E-02,  8.065795724472015E-04);
        vec3 neptun_vel(7.929670145395670E-04, 3.049616262068798E-03, -8.104050065051757E-05);
        vec3 saturn_vel(5.202769599394199e-3,8.527763066533207e-4,-2.216361381984069e-4);
        vec3 mars_vel(5.680349952465525E-04, 1.518932715814366E-02, 3.043415082859981E-04);
        vec3 jupiter_vel(6.472894738201937e-03, -3.365002410907006e-03, -1.308549331239642e-04);
        //vec3 earth_vel(-6.828707544910924e-03, 1.585182505936633e-02, 8.139769320198632e-08);
        //vec3 mercury_vel(2.218153157209631E-02, 5.607002210997516E-03, -1.576726276143236E-03);

        vec3 mercury_vel(0,12.44,0);
        vec3 earth_vel(0, 2*M_PI, 0);




        pluto_vel = pluto_vel*365.25;
        uranus_vel = uranus_vel*365.25;
        venus_vel = venus_vel*365.25;
        neptun_vel = neptun_vel*365.25;
        saturn_vel = saturn_vel*365.25;
        mars_vel = mars_vel*365.25;
        //earth_vel = earth_vel*365.25;
        jupiter_vel = jupiter_vel*365.25;
        //mercury_vel = mercury_vel*365.25;




        int N = 1e4 ;


        //initialzing bodies
        body Pluto(6.545e-9, pluto_pos, pluto_vel, "pluto");
        body Venus(0.0000024335, venus_pos, venus_vel, "venus");
        body Uranus(0.000043405, uranus_pos, uranus_vel, "uranus");
        body Neptun(0.0000512, neptun_pos, neptun_vel, "neptun");
        body Mars(3.213e-7, mars_pos, mars_vel,"mars");
        body Saturn(0.0002857, saturn_pos, saturn_vel,"saturn");
        body Earth(0.000003, earth_pos, earth_vel, "earth");
        body Jupiter(0.000949, jupiter_pos,jupiter_vel, "jupiter");
        body Mercury(1.6427e-7, mercury_pos,mercury_vel, "Mercury");


        double final_time = 10;
        double dt = final_time/((double) N);
        quantities solar_system;

        //adding bodies
        solar_system.add_new_body(Pluto, store_values(N));
        solar_system.add_new_body(Venus, store_values(N));
        solar_system.add_new_body(Uranus, store_values(N));
        solar_system.add_new_body(Mars, store_values(N));
        solar_system.add_new_body(Neptun, store_values(N));
        solar_system.add_new_body(Saturn, store_values(N));
        solar_system.add_new_body(Earth, store_values(N));
        solar_system.add_new_body(Jupiter, store_values(N));
        solar_system.add_new_body(Mercury, store_values(N));


        vec3 sun_vel = {0,0,0};
        vec3 sun_pos = {0,0,0};


        //this ensures that the total momentum is equal zero in the begining
        for(int v = 0; v < solar_system.list_of_bodies.size(); v++){
            body &current_planet = solar_system.list_of_bodies[v];
            sun_vel += current_planet.velocity*current_planet.mass*(-1);
            sun_pos += current_planet.position*current_planet.mass;
        }


        vec3 P;
        double eps = 1e-8;
        //extra check for momentum
        for(int planets = 0; planets < solar_system.list_of_bodies.size(); planets++){
            body &this_planet = solar_system.list_of_bodies[planets];
            P += this_planet.mass*this_planet.velocity;
        }
        if((P + sun_vel).length() >= eps ){
            cout << P + sun_vel << endl;
            cout << "total momentum is not conserved in the begining, please try agian" << endl;
            exit(1);


            body Sun(1.0, sun_pos, sun_vel, "sun");
            solar_system.add_new_body(Sun, store_values(N));








            solver solutions(dt, solar_system);

            double beta = 2.0;

            vec3 p_start;
            vec3 p_end;

            clock_t start_time, end_time;
            double time_used;
            start_time = clock();

            p_start = solar_system.list_of_bodies[0].position;
            //main loop
            for(int i=0; i<N; i++){
                solutions.forward_euler(solar_system, i, beta);

            }


            p_end = solar_system.list_of_bodies[0].position;

            cout << (fabs(p_end.length()-p_start.length())/(p_start.length())) << endl;
            end_time = clock();
            solar_system.write_to_file(beta);

            time_used = ((double)(end_time - start_time))/CLOCKS_PER_SEC;
            cout << "time used for N = " << N << "and dt = " << dt << "is t = " << time_used << endl;
            //cout << store_values.theta[0] - store_values.theta[int(N/100)] << endl;

            return 0;

        }

}
