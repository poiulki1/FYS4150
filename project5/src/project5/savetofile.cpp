#include "savetofile.h"

savetofile::savetofile(string fname, string number_arg):
    filename(fname),
    number_arg(number_arg)
{
    outfile.open("variance.txt");

}

void savetofile::save(string arg, int b, mat &M, rowvec &y, int run_length){
    ofstream file;

    if(arg == "mc"){
        file.open("../result/" + filename + "_mc_b_" + to_string(b) +"_acddief_"+ number_arg+ ".txt");
    }
    else if(arg == "rk4"){
        file.open("../result/" + filename + "_rk4_b_" + to_string(b) + "_acddief_" + number_arg + ".txt");
    }
    //else{
    //    file.open("../result/" + "variance_b" + to_string(b) + "_acddief_" + number_arg + ".txt");
    //}


    for(int i = 0; i < run_length; i++){
        file <<setw(13) <<y(i) << setw(13) <<  M(i, 0) <<
               setw(13) <<  M(i, 1) <<
               setw(13) <<  M(i, 2) << endl;
    }


    if(arg == "mc"){
        //file.open(filename + "_mc_b_" + to_string(b) +"_acddief_"+ number_arg+ ".txt");
         cout << filename +"_mc_b" + to_string(b) + "_acddief_" + number_arg + ".txt is written" << endl;
         cout << "----------------------------------------------------------------- "<< endl;
    }
    else if(arg == "rk4"){
        cout << filename + "_rk4_b_" + to_string(b) + "_acddief_" + number_arg + ".txt is written" << endl;
        cout << "----------------------------------------------------------------- "<< endl;
    }

    file.close();



}

void savetofile::save_each_step(double M1, double M2, double M3, double y){

    outfile <<setw(13) <<  y  <<
              setw(13) <<  M1 <<
              setw(13) <<  M2 <<
              setw(13) <<  M3 << endl;

}
