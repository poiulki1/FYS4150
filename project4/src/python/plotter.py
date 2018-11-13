import numpy as np
import matplotlib.pyplot as plt

def analytic_cv(temp):
    k = 1.0
    J = 1.0
    beta = 1.0/(temp*k)
    beta2 = 1.0/(temp*temp*k)
    e_arg = 8*beta*J
    fac1 = 128.0*np.exp(-e_arg)+ 128.0*np.exp(e_arg)
    fac2 = 16*J*np.exp(-e_arg) - 16*J*np.exp(e_arg)
    dem = 2*np.exp(e_arg) + 2*np.exp(-e_arg)+12.0
    cv = beta2*((fac1/dem)-((fac2/dem)**2))*(1/4.0)    

    return cv
    
def plot_4b_cv(filename):
    data = np.loadtxt(filename)

    temp = data[:40,1]
    cv_analytic = analytic_cv(temp)
    cv4 = data[:40,5]
    cv5 = data[40:80,5]
    cv6 = data[80:120,5]
    cv7 = data[120:160,5]
    plt.figure()
    plt.title("Heat capacity per spin with $10^{4}$ MC cycles")
    plt.xlabel("Temperature [K]")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.plot(temp, cv_analytic,"--", label="analytic")
    plt.plot(temp, cv4, label="numeric")
    plt.legend(loc="best")
    plt.savefig("4b_cv_mc4.pdf")
    plt.figure()
    plt.title("Heat capacity per spin with $10^{5}$ MC cycles")
    plt.xlabel("Temperature [K]")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.plot(temp, cv_analytic,"--", label="analytic")
    plt.plot(temp, cv5, label="numeric")
    plt.legend(loc="best")
    plt.savefig("4b_cv_mc5.pdf")
    plt.figure()
    plt.title("Heat capacity per spin with $10^{6}$ MC cycles")
    plt.xlabel("Temperature [K]")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.plot(temp, cv_analytic, "--",label="analytic")
    plt.plot(temp, cv6, label="numeric")
    plt.legend(loc="best")
    plt.savefig("4b_cv_mc6.pdf")
    plt.figure()
    plt.title("Heat capacity per spin with $10^{7}$ MC cycles")
    plt.xlabel("Temperature [K]")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.plot(temp, cv_analytic,"--", label="analytic")
    plt.plot(temp, cv7, label="numeric")
    plt.legend(loc="best")
    plt.savefig("4b_cv_mc7.pdf")
    
def plot_energy_dist(temp1, temp24):
    energy_dist_temp1 = np.loadtxt(temp1)
    energy_dist_temp24 = np.loadtxt(temp24)

    energy_temp1, cnt_temp1 = np.unique(energy_dist_temp1, return_counts=True)
    energy_temp24, cnt_temp24 = np.unique(energy_dist_temp24, return_counts=True)

    plt.figure()
    plt.title("Energy distribution, $20$x$20$ lattice, total MC cycles = $10^{6}$ ")
    plt.xlabel("Energy [J]")
    plt.ylabel("P(E)")
    A1, B1 = np.histogram(energy_dist_temp1)
    h1 =  A1/float(sum(A1))
    width1 = B1[1]-B1[0]
    plt.bar(B1[:len(B1)-1],h1, width=width1)
    plt.savefig("4D_energy_dist1.pdf")
    
    plt.figure()
    plt.title("Energy distribution, $20$x$20$ lattice, total MC cycles = $10^{6}$ ")
    plt.xlabel("Energy [J]")
    plt.ylabel("P(E)")
    A24, B24 = np.histogram(energy_dist_temp24)
    h24 =  A24/float(sum(A24))
    width24 = B24[1]-B24[0]
    plt.bar(B24[:len(B24)-1],h24, width=width24)
    plt.savefig("4D_energy_dist24.pdf")
    #plt.show()

def plot_accepted_temp(accepted_random, accepted_not_random):
    random = np.loadtxt(accepted_random)
    not_random = np.loadtxt(accepted_not_random)
    
    temp = random[:,1]
    
    acc_random = random[:,0]
    acc_not_random = not_random[:,0]

    plt.figure()
    plt.xlabel("Temperature [K]")
    plt.ylabel("Accepted spin flips/total mc cycles")
    plt.plot(temp,acc_random, label="Disordered initial state")
    plt.plot(temp, acc_not_random, label="Ordered initial state")
    plt.title("Accepted spin flips per MC cycle per spin")
    plt.legend(loc="best")
    plt.savefig("4C_accepted_temp.pdf")

def plot_accepted_mc(filename_random1, filename_random24, filename_not_random1, \
                     filename_not_random24):
    
    acc_random_temp1 = np.loadtxt(filename_random1)
    acc_random_temp24 = np.loadtxt(filename_random24)
    
    acc_not_random_temp1 = np.loadtxt(filename_not_random1)
    acc_not_random_temp24 = np.loadtxt(filename_not_random24)

    mc_steps = acc_random_temp1[:,0]

    plt.figure()
    plt.title("Accepted spin flips")
    plt.xlabel("MC step")
    plt.plot(mc_steps, acc_random_temp1[:,1], label="Disordered initial state")
    plt.plot(mc_steps, acc_not_random_temp1[:,1], label="Ordered initial state")
    plt.legend()
    plt.savefig("4C_accepted_mc_temp1.pdf")
    plt.figure()
    plt.title("Accepted spin flips")
    plt.xlabel("MC step")
    plt.plot(mc_steps, acc_random_temp24[:,1], label="Disordered initial state")
    plt.plot(mc_steps, acc_not_random_temp24[:,1], label="Ordered initial state")
    plt.legend()
    plt.savefig("4C_accepted_mc_temp24.pdf")
    
    
def plot_mean(filename1, filename2, filename3, filename4):
    mean_temp1_random = np.loadtxt(filename1)
    mean_temp24_random = np.loadtxt(filename2)

    mean_temp1_not_random = np.loadtxt(filename3)
    mean_temp24_not_random = np.loadtxt(filename4)

    mc_steps = mean_temp1_random[:,0]

    energy_temp1_random = mean_temp1_random [:,1]
    mag_temp1_random = mean_temp1_random[:,2]
    energy_temp24_random = mean_temp24_random[:,1]
    mag_temp24_random = mean_temp24_random[:,2]

    energy_temp1_not_random = mean_temp1_not_random [:,1]
    mag_temp1_not_random = mean_temp1_not_random[:,2]
    energy_temp24_not_random = mean_temp24_not_random[:,1]
    mag_temp24_not_random = mean_temp24_not_random[:,2]

    plt.figure()
    plt.plot(mc_steps, energy_temp1_random)
    plt.title("Mean energy (Disordered initial state)")
    plt.xlabel("MC step")
    plt.ylabel("Energy [J]")
    plt.savefig("mean_E_random_temp1.pdf")

    plt.figure()
    plt.plot(mc_steps, mag_temp1_random)
    plt.title("Mean magnetisation (Disordered initial state)")
    plt.xlabel("MC step")
    plt.ylabel("Magnetisation [J/T]")
    plt.savefig("mean_M_random_temp1.pdf")

    plt.figure()
    plt.plot(mc_steps, energy_temp24_random)
    plt.title("Mean energy (Disordered initial state)")
    plt.xlabel("MC step")
    plt.ylabel("Energy [J]")
    plt.savefig("mean_E_random_temp24.pdf")

    plt.figure()
    plt.plot(mc_steps, mag_temp24_random)
    plt.title("Mean magnetisation (Disordered initial state)")
    plt.xlabel("MC step")
    plt.ylabel("Magnetisation [J/T]")
    plt.savefig("mean_M_random_temp24.pdf")
    
    plt.figure()
    plt.plot(mc_steps, energy_temp1_not_random)
    plt.title("Mean energy (Ordered initial state)")
    plt.xlabel("MC step")
    plt.ylabel("Energy [J]")
    plt.savefig("mean_E_not_random_temp1.pdf")
    
    plt.figure()
    plt.plot(mc_steps, mag_temp1_not_random)
    plt.title("Mean magnetisation (Ordered initial state)")
    plt.xlabel("MC step")
    plt.ylabel("Magnetisation [J/T]")
    plt.savefig("mean_M_not_random_temp1.pdf")
    
    plt.figure()
    plt.plot(mc_steps, energy_temp24_not_random)
    plt.title("Mean energy (Ordered initial state)")
    plt.xlabel("MC step")
    plt.ylabel("Energy [J]")
    plt.savefig("mean_E_not_random_temp24.pdf")
    
    plt.figure()
    plt.plot(mc_steps, mag_temp24_not_random)
    plt.title("Mean magnetisation (Ordered initial state)")
    plt.xlabel("MC step")
    plt.ylabel("Magnetisation [J/T]")
    plt.savefig("mean_M_not_random_temp24.pdf")



#calling
#plot_4b_cv("4B_data.txt")
#plot_energy_dist("energy_dist_temp1.txt", "energy_dist_temp24.txt")

#plot_accepted_temp("accepted_temp_random.txt", "accepted_temp_not_random.txt")

#plot_accepted_mc("accepted_mc_random_temp1.txt", \
#                 "accepted_mc_random_temp24.txt", \
#                 "accepted_mc_not_random_temp1.txt", \
#                 "accepted_mc_not_random_temp24.txt")

plot_mean("4C_temp1_random.txt", "4C_temp24_random.txt", \
                 "4C_temp1_not_random.txt", "4C_temp24_not_random.txt")
