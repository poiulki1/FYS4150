import numpy as np
import matplotlib.pyplot as plt


argument1 = "earth"
argument2 = "earthverlet_beta_"
info = "_entire_quantity_info_verlet"
exstension = ".txt"
bodies = ["sun", "mercury", "venus", "earth", "mars",
"jupiter", "saturn", "uranus", "neptun", "pluto"]
number_of_files = int(input("number of files: "))
number_of_years = int(input("number of years: "))
beta = [2.0,2.3,2.5, 2.8,3.0]
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

data = []
#for i in range(0,number_of_files):
#    data.append(np.loadtxt(bodies[i] +  exstension))
    #p2 = np.loadtxt(bodies[i] + str(i) + exstension)



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
    plt.ylabel("Energy [$\\frac{AU^{2}}{M_{\\odot}yr^{2}}}$]")
    plt.grid("on")
    #plt.savefig("energy_in_solar_system.pdf")
    plt.show()

    plt.plot(time, total_angular_momentum[0]- total_angular_momentum)

    plt.title("Total total angular momentum for the solar system")
    #plt.legend(["Total energy", " Kinetic energy", "Potential energy"])
    plt.xlabel("time [yrs]")
    plt.ylabel("angular momentum[$\\frac{Au^{2}}{yrM_{\\odot}}}$]")
    plt.grid("on")
    #plt.savefig("angular_momentum_solar_system.pdf")
    plt.show()


    plt.plot(time, total_momentum[0]-total_momentum)
    plt.title("Total momnetum for the entire solar system")
    #plt.legend(["Total energy", " Kinetic energy", "Potential energy"])
    plt.xlabel("time [yrs]")
    plt.ylabel("momentum[$\\frac{AU}{M_{\\odot}yr}$]")
    plt.grid("on")
    #plt.savefig("momentum_solar_system.pdf")
    plt.show()




#solar_system_info(bodies, info, exstension, number_of_files, number_of_years)
def earth_sun_num_comparison(data1,data2, number_of_files):


    for j in range(number_of_files):
        plt.plot(data1[j][:,0], data1[j][:,1])
        #print(data1[j][:,0])
        #plt.plot(data2[j][:,0], data2[j][:,1])


    #plt.plot(0,0, "-ok")

    plt.ylabel("Energy")
    plt.xlabel("$time[yr]$")
    plt.legend(["total kinetic energy",
    "total potential energy", "total energy"], loc = "upper right")
    plt.title("Different energy quantities, Forward Euler")
    plt.savefig("solar_system_almost_enitre.pdf")
    plt.show()



#earth_sun_num_comparison(data, 2, number_of_files)
#energy_plot(number_of_files,bodies, kinetic_and_potential, exstension)



def energy_earth_sun():
    s_vv = np.loadtxt("sun_quantity_info_verlet_final.txt")
    e_vv = np.loadtxt("earth_quantity_info_verlet_final.txt")

    s_fe = np.loadtxt("sun_quantity_info_euler_final.txt")
    e_fe = np.loadtxt("earth_quantity_info_euler_final.txt")

    mv1 = s_fe[:,0] + e_fe[:,2]
    #K1 = abs(K1 - K1[0])/abs(K1[0])

    L1 = s_fe[:,1] + e_fe[:,3]
    #P1 = abs(P1 - P1[0])/abs(P1[0])

    mv2 = s_vv[:,0] + e_vv[:,2]
    #K2 = abs(K2 - K2[0])/abs(K2[0])

    L2 = s_vv[:,1] + e_vv[:,3]
    #P2 = abs(P2 - P2[0])/abs(P2[0])

    #total_energy1 = K1 + P1
    #total_energy2 = K2 + P2

    time = np.linspace(0,10,len(mv1))

    plt.title("The earth sun angular and linear momentum, Forward Euler")

    plt.plot(time,mv1, "-r")
    plt.plot(time,L1, "-k")
    plt.xlabel("Time [yrs]")
    plt.ylabel("Energy")
    #plt.plot(time, (total_energy1[0] - total_energy1),  "-b")
    plt.legend(["$|\\mathbf{P}| = m|\\mathbf{v}|$", "$|\\mathbf{L}| = |\\mathbf{r} x \\mathbf{v}|$"])
    #plt.axis([xmin,xmax,ymin,ymax])
    plt.savefig("euler_sun_earth_momentum2.pdf")


    plt.show()


    plt.title("The earth sun angular and linear momentum, Velocity Verlet")
    plt.xlabel("Time [yrs]")
    plt.ylabel("Energy")
    plt.plot(time,mv2, "-r")
    plt.plot(time,L2, "-k")

    #plt.plot(time, total_energy2[0] - total_energy2,  "-b")
    plt.legend(["$|\\mathbf{P}| = m|\\mathbf{v}|$","$|\\mathbf{L}| = |\\mathbf{r} x \\mathbf{v}|$"])
    plt.savefig("verlet_sun_earth_momentum.pdf")
    plt.show()

    #plt.ylabel("Energy")
    #plt.xlabel("$time[yr]$")
    #plt.title("The total energy in the system, simulated for 100yrs")
    #plt.legend(["kinetic energy FE",
    #"kinetic energy VV"
    #"potential energy FE", "potential energy VV"], loc = "upper right")

    #plt.savefig("earth_sun_total_energy_verlet.pdf")
    #plt.show()
#energy_earth_sun()
def perhelion_plot():
    data = np.loadtxt("mercuryperihelion_new.txt")
    data2 = np.loadtxt("sunperihelion_new.txt")
    #print(data)

    angles = []

    x = data[:,0]
    y= data[:,1]

    x0 = data2[:,0]
    y0 = data2[:,1]

    for i in range(2, len(x)):
        prev_prev_R = np.sqrt((x[i-2] - x0[i-2])**2 + (y[i-2] - y0[i-2])**2)
        prev_R = np.sqrt((x[i-1] - x0[i-1])**2 + (y[i-1] - y0[i-1])**2)
        R = np.sqrt((x[i] - x0[i])**2 + (y[i] - y0[i])**2)

        if R > prev_R and prev_R < prev_prev_R:
            angles.append(np.arctan2(y[i-1] - y0[i-1], x[i-1] - x0[i-1]))


    angles = np.array(angles)
    print(angles - angles[0])
    plt.plot(np.linspace(0,100,len(angles)), (angles - angles[0]))
    plt.title("The perihelion angle measured during the course of 100 yrs")
    plt.xlabel("Time[yrs]")
    plt.ylabel("$\\theta_{p}$[Radians]")
    #plt.savefig("perihelion_angle_plot.pdf")
    #plt.show()

perhelion_plot()

def earth_sun_jup():
    earth1 = np.loadtxt("earthsun_earth_jup_fixed_com1000.txt")
    jupiter1 = np.loadtxt("jupitersun_earth_jup_fixed_com1000.txt")
    sun1 = np.loadtxt("sunsun_earth_jup_fixed_com1000.txt")


    plt.plot(earth1[:,0], earth1[:,1])
    plt.plot(jupiter1[:,0], jupiter1[:,1])
    plt.plot(sun1[:,0], sun1[:,1])
    plt.xlabel("x[AU]")
    plt.ylabel("y[AU]")
    plt.title("Earth, Sun, Jupiter system with 1000 times jupiter mass")
    plt.legend(["earth", "jupiter", "sun"])
    plt.savefig("earth_sun_jupiter_m1000.pdf")
    plt.show()

#earth_sun_jup()
#perhelion_plot()












#earth_sun_num_comparison(np.loadtxt("mercury.txt"),2 ,number_of_files)
