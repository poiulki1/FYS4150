import numpy as np
import matplotlib.pyplot as plt


argument1 = "earth"
argument2 = "earthverlet_beta_"
info = "_quantity_info"
exstension = "_perhelion.txt"
bodies = ["sun", "mercury", "earth", "jupiter", "mars", "saturn"]
number_of_files = int(input("number of files: "))
number_of_years = int(input("number of years: "))
beta = [2.0,2.3,2.5,2.8,3.0]
euler = []
verlet = []


"""m = np.loadtxt("system_pot_energy.txt")
n = np.loadtxt("system_kin_energy.txt")

#print(m+n)
#print(np.shape(m), np.shape(n))

plt.plot(np.linspace(0,30,len(m)), m)
plt.plot(np.linspace(0,30,len(m)), n)
plt.show()
exit()"""

"""data = []
for i in range(0,number_of_files):
    data.append(np.loadtxt(bodies[i] + exstension))
    #p2 = np.loadtxt(bodies[i] + str(i) + exstension)

"""

def time_plot():
    dt = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
    
    euler_time = [9e-5, 5.05e-4, 7.91e-3, 5.24e-2, 5.33e-1, 5.17,  58]
    vv_time = [1.11e-4, 1.07e-3, 8.93e-3, 9.71e-2, 7.86e-1, 7.73, 88.4]

    euler_error = [4.51e1, 1.43, 4.92e-1, 7.32e-2, 7.8e-3, 7.82e-4, 7.83e-5]
    vv_error = [1.85e-1, 1.53e-5, 6.78e-8, 8.73e-8, 8.88e-8, 8.88e-8,8.88e-8] 
    
    plt.plot(np.log(dt), np.log(euler_error), label="Forward Euler")
    plt.plot(np.log(dt), np.log(vv_error), label="Velocity Verlet")
    plt.xlabel("$\Delta t$")
    plt.ylabel("Relative error in position")
    plt.title("Precision comparison")
    plt.legend()
    plt.show()
    
def jupiter_plot():
    jup0 = np.loadtxt("jupiter0.txt")
    earth0 = np.loadtxt("earth0.txt")
    sun0 = np.loadtxt("sun0.txt")
    
    jup10 = np.loadtxt("jupiter10.txt")
    earth10 = np.loadtxt("earth10.txt")
    sun10 = np.loadtxt("sun10.txt")

    jup1000 = np.loadtxt("jupiter1000.txt")
    earth1000 = np.loadtxt("earth1000.txt")
    sun1000 = np.loadtxt("sun1000.txt")

    plt.figure()
    plt.plot(jup0[:,0], jup0[:,1], label="Jupiter")
    plt.plot(sun0[:,0], sun0[:,1], label="Sun")
    plt.plot(earth0[:,0], earth0[:,1], label="Earth")
    plt.xlabel("X [AU]")
    plt.ylabel("Y [AU]")
    plt.legend()
    plt.title("Sun-Jupiter-Earth system with real Jupiter mass")
    plt.savefig("jupiter_0.pdf")

    plt.figure()
    plt.plot(jup10[:,0], jup10[:,1], label="Jupiter")
    plt.plot(sun10[:,0], sun10[:,1], label="Sun")
    plt.plot(earth10[:,0], earth10[:,1], label="Earth")
    plt.xlabel("X [AU]")
    plt.ylabel("Y [AU]")
    plt.legend()
    plt.title("Sun-Jupiter-Earth system with Jupiter mass increased by factor of 10")
    plt.savefig("jupiter_10.pdf")
    
    plt.figure()
    plt.plot(jup1000[:,0], jup1000[:,1], label="Jupiter")
    plt.plot(sun1000[:,0], sun1000[:,1], label="Sun")
    plt.plot(earth1000[:,0], earth1000[:,1], label="Earth")
    plt.xlabel("X [AU]")
    plt.ylabel("Y [AU]")
    plt.axis([-5, 15, -15, 5])
    plt.title("Sun-Jupiter-Earth system with Jupiter mass increased by factor of 1000")
    plt.legend()
    plt.savefig("jupiter_1000.pdf")

def solar_system_info(bodies, info, exstension, number_of_files,t):
    quantities = []

    for i in range(number_of_files):
        quantities.append(np.loadtxt(bodies[i] + info + exstension))

    total_kinetic_energy = 0
    total_potential_energy = 0
    total_angular_momentum = 0
    total_momentum = 0

    for j in range(number_of_files):
        total_kinetic_energy += quantities[j][:,0]
        total_potential_energy += quantities[j][:,1]
        total_angular_momentum += quantities[j][:,3]
        total_momentum += quantities[j][:,2]

    time = np.linspace(0,t, len(total_kinetic_energy))

    plt.plot(time, total_kinetic_energy + total_potential_energy)
    plt.plot(time, total_kinetic_energy)
    plt.plot(time, total_potential_energy)
    plt.title("Total energy for the entire solar system")
    plt.legend(["Total energy", " Kinetic energy", "Potential energy"])
    plt.xlabel("time [yrs]")
    plt.ylabel("Energy [$J$]")
    plt.grid("on")
    #plt.savefig("energy_in_solar_system.pdf")
    plt.show()

    plt.plot(time, abs(total_angular_momentum[0]-total_angular_momentum)/abs(total_angular_momentum[0]))

    plt.title("Total total angular momentum for the solar system")
    #plt.legend(["Total energy", " Kinetic energy", "Potential energy"])
    plt.xlabel("time [yrs]")
    plt.ylabel("angular momentum")
    plt.grid("on")
    plt.savefig("angular_momentum_solar_system.pdf")
    plt.show()


    plt.plot(time, abs(total_momentum[0]-total_momentum)/abs(total_momentum[0]))
    plt.title("Total momnetum for the entire solar system")
    #plt.legend(["Total energy", " Kinetic energy", "Potential energy"])
    plt.xlabel("time [yrs]")
    plt.ylabel("momentum]")
    plt.grid("on")
    plt.savefig("momentum_solar_system.pdf")
    plt.show()




#solar_system_info(number_of_files, number_of_years)
def earth_sun_num_comparison(data1,data2, number_of_files):


    #for j in range(number_of_files):
    #    plt.plot(data1[j][:,0], data1[j][:,1], "-r")
        #plt.plot(data2[j][:,0], data2[j][:,1])


    #plt.plot(0,0, "-ok")
    plt.plot(np.linspace(0,100,1e7), data1)
    plt.ylabel("$\\theta_{p}$")
    plt.xlabel("$time$[yrs]")
    plt.legend(["2"], loc = "upper right")
    plt.title("The perhelion angle for mercury-sun system (100 earth yrs with $N = 1e7$)")
    #plt.savefig("merc_sun_perhelion_angle_100yrs_einstein.pdf")
    plt.show()



#energy_plot(number_of_files,bodies, kinetic_and_potential, exstension)


def perhelion_plot():
    data = np.loadtxt("mercury.txt")
    x = data[:,0]
    y = data[:,1]

    angles = []
    for i in range(2, len(x)):
        prev_prev_R = np.sqrt(x[i-2]*x[i-2] + y[i-2]*y[i-2])
        prev_R = np.sqrt(x[i-1]*x[i-1] + y[i-1]*y[i-1])
        R = np.sqrt(x[i]*x[i] + y[i]*y[i])

        if R > prev_R and prev_R < prev_prev_R:
            angles.append(np.arctan(y[i-1]/x[i-1]))

    angles = np.array(angles)
    plt.plot(np.linspace(0,100,len(angles)), angles)
    plt.show()


#jupiter_plot()
#perhelion_plot()


#earth_sun_num_comparison(np.loadtxt("mercury.txt"),2 ,number_of_files)
