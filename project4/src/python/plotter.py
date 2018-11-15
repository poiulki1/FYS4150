import numpy as np
import matplotlib.pyplot as plt

def max_point(arr):
    #temp = np.loadtxt("4E_data_L40_zoom.txt")[:,1]
    temp = np.loadtxt("4E_data_L40.txt")[:,1]
    x_max = temp[np.argmax(arr)]
    y_max = arr.max()
    return x_max, y_max

def plot_4e(arg="not_zoom"):
    if(arg=="zoom"):
        filename40 = "4E_data_L40_zoom.txt"
        filename60 = "4E_data_L60_zoom.txt"
        filename80 = "4E_data_L80_zoom.txt"
        filename100 = "4E_data_L100_zoom.txt"
    else:       
        filename40 = "4E_data_L40.txt"
        filename60 = "4E_data_L60.txt"
        filename80 = "4E_data_L80.txt"
        filename100 = "4E_data_L100.txt"

    data40 = np.loadtxt(filename40)
    data60 = np.loadtxt(filename60)
    data80 = np.loadtxt(filename80)
    data100 = np.loadtxt(filename100)
    
    temp = data40[:,1]

    mean_E40 = data40[:,2]
    mean_E60 = data60[:,2]
    mean_E80 = data80[:,2]
    mean_E100 = data100[:,2]
    
    mean_Mabs40 = data40[:,8]
    mean_Mabs60 = data60[:,8]
    mean_Mabs80 = data80[:,8]
    mean_Mabs100 = data100[:,8]
    
    cv40 = data40[:,5]
    cv60 = data60[:,5]
    cv80 = data80[:,5]
    cv100 = data100[:,5]

    xi40 = data40[:,10]
    xi60 = data60[:,10]
    xi80 = data80[:,10]
    xi100 = data100[:,10]

    plt.figure()
    plt.title("Mean energy per spin")
    plt.plot(temp, mean_E40, label="L=40")
    plt.plot(temp, mean_E60, label="L=60")
    plt.plot(temp, mean_E80, label="L=80")
    plt.plot(temp, mean_E100, label="L=100")
    plt.ylabel("<E> [J]")
    plt.xlabel("Temperature [K]")
    plt.legend(loc="best")
    
    plt.figure()
    plt.title("Mean absolute magnetisation per spin")
    plt.plot(temp, mean_Mabs40, label="L=40")
    plt.plot(temp, mean_Mabs60, label="L=60")
    plt.plot(temp, mean_Mabs80, label="L=80")
    plt.plot(temp, mean_Mabs100, label="L=100")
    plt.ylabel("<|M|> [J/T]")
    plt.xlabel("Temperature [K]")
    plt.legend(loc="best")
    
    plt.figure()
    plt.title("Heat capacity per spin")
    plt.plot(temp, cv40, label="L=40")
    plt.plot(max_point(cv40)[0], max_point(cv40)[1], "o")
    plt.plot(temp, cv60, label="L=60")
    plt.plot(max_point(cv60)[0], max_point(cv60)[1], "o")
    plt.plot(temp, cv80, label="L=80")
    plt.plot(max_point(cv80)[0], max_point(cv80)[1], "o")
    plt.plot(temp, cv100, label="L=100")
    plt.plot(max_point(cv100)[0], max_point(cv100)[1], "o")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.xlabel("Temperature [K]")
    plt.legend(loc="best")
    
    plt.figure()
    plt.title("Magnetic susceptibility per spin")
    plt.plot(temp, xi40, label="L=40")
    plt.plot(max_point(xi40)[0], max_point(xi40)[1], "o")
    plt.plot(temp, xi60, label="L=60")
    plt.plot(max_point(xi60)[0], max_point(xi60)[1], "o")
    plt.plot(temp, xi80, label="L=80")
    plt.plot(max_point(xi80)[0], max_point(xi80)[1], "o")
    plt.plot(temp, xi100, label="L=100")
    plt.plot(max_point(xi100)[0], max_point(xi100)[1], "o")
    plt.ylabel("$\\chi$ [$Js^{2}$]")
    plt.xlabel("Temperature [K]")
    plt.legend(loc="best")
    plt.show()

def analytic_cv(temp):
    k = 1.0
    J = 1.0
    beta = 1.0/(temp*k)
    beta2 = 1.0/(temp*temp*k)
    e_arg = 8*beta*J
    fac1 = 128.0*np.exp(-e_arg)+ 128.0*np.exp(e_arg)
    fac2 = 16*J*np.exp(-e_arg) - 16*J*np.exp(e_arg)
    dem = 2*np.exp(e_arg) + 2*np.exp(-e_arg)+12.0
    cv = beta2*((fac1/dem)-((fac2/dem)**2))*(float(1./4.0))    

    return cv
    
def plot_4b_cv(filename):
    data = np.loadtxt(filename)

    temp = data[:40,1]
    cv_analytic = analytic_cv(temp)
    
    cv4 = data[:40,5]
    cv5 = data[40:80,5]
    cv6 = data[80:120,5]
    #cv7 = data[120:160,5]
    plt.figure()
    plt.title("Heat capacity per spin with $10^{4}$ MC cycles")
    plt.xlabel("Temperature [K]")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.plot(temp, cv_analytic,"--", label="analytic")
    plt.plot(temp, cv4, label="numeric")
    plt.legend(loc="best")
    #plt.savefig("4b_cv_mc4.pdf")
    plt.figure()
    plt.title("Heat capacity per spin with $10^{5}$ MC cycles")
    plt.xlabel("Temperature [K]")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.plot(temp, cv_analytic,"--", label="analytic")
    plt.plot(temp, cv5, label="numeric")
    plt.legend(loc="best")
    #plt.savefig("4b_cv_mc5.pdf")
    plt.figure()
    plt.title("Heat capacity per spin with $10^{6}$ MC cycles")
    plt.xlabel("Temperature [K]")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.plot(temp, cv_analytic, "--",label="analytic")
    plt.plot(temp, cv6, label="numeric")
    plt.legend(loc="best")
    #plt.savefig("4_cv_mc6.pdf")
    """
    plt.figure()
    plt.title("Heat capacity per spin with $10^{7}$ MC cycles")
    plt.xlabel("Temperature [K]")
    plt.ylabel("$C_{v}$ [J/K]")
    plt.plot(temp, cv_analytic,"--", label="analytic")
    plt.plot(temp, cv7, label="numeric")
    plt.legend(loc="best")"""
    #plt.savefig("4b_cv_mc7.pdf")
    plt.show()
    
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
    #plt.savefig("4D_energy_dist1.pdf")
    
    plt.figure()
    plt.title("Energy distribution, $20$x$20$ lattice, total MC cycles = $10^{6}$ ")
    plt.xlabel("Energy [J]")
    plt.ylabel("P(E)")
    A24, B24 = np.histogram(energy_dist_temp24)
    h24 =  A24/float(sum(A24))
    width24 = B24[1]-B24[0]
    plt.bar(B24[:len(B24)-1],h24, width=width24)
    #plt.savefig("4D_energy_dist24.pdf")
    plt.show()

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
plot_4e()#arg="zoom")

#plot_4b_cv("TEST_4b_MPI.txt")
#plot_4b_cv("4B_data.txt")
#plot_energy_dist("energy_dist_temp1_MPI.txt", "energy_dist_temp24_MPI.txt")

#plot_accepted_temp("accepted_temp_random.txt", "accepted_temp_not_random.txt")

#plot_accepted_mc("accepted_mc_random_temp1.txt", \
#                 "accepted_mc_random_temp24.txt", \
#                 "accepted_mc_not_random_temp1.txt", \
#                 "accepted_mc_not_random_temp24.txt")

#plot_mean("4C_temp1_random.txt", "4C_temp24_random.txt", \
#                 "4C_temp1_not_random.txt", "4C_temp24_not_random.txt")
