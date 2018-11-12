import numpy as np
import matplotlib.pyplot as plt


def plot_accepted_temp(accepted_random, accepted_not_random):
    random = np.loadtxt(accepted_random)
    not_random = np.loadtxt(accepted_not_random)
    
    temp = random[:,1]
    
    acc_random = random[:,0]
    acc_not_random = not_random[:,0]

    plt.figure()
    plt.plot(temp,acc_random, label="Disordered initial state")
    plt.plot(temp, acc_not_random, label="Ordered initial state")
    plt.title("Accepted spin flips")
    plt.legend()

    plt.figure()
    plt.plot(mc_steps, accepted_random_temp1)
    plt.title()

    plt.figure()
    plt.plot(mc_steps, accepted_random_temp24)
    plt.show()

def plot_accepted_mc(filename_random1, filename_random24, filename_not_random1, \
                     filename_not_random24):
    
    acc_random_temp1 = np.loadtxt(filename_random1)
    acc_random_temp24 = np.loadtxt(filename_random24)
    
    acc_not_random_temp1 = np.loadtxt(filename_not_random1)
    acc_not_random_temp24 = np.loadtxt(filename_not_random24)

    mc_steps = acc_random_temp1[:,0]

    plt.figure()
    plt.title("Accepted spins, temperature=1 [kT/J]")
    plt.plot(mc_steps, acc_random_temp1[:,1], label="Disordered initial state")
    plt.plot(mc_steps, acc_not_random_temp1[:,1], label="Ordered initial state")
    plt.legend()
    
    plt.figure()
    plt.title("Accepted spins, temperature=2.4 [kT/J])")
    plt.plot(mc_steps, acc_random_temp24[:,1], label="Disordered initial state")
    plt.plot(mc_steps, acc_not_random_temp24[:,1], label="Ordered initial state")
    plt.legend()
    
    plt.show()
    
    
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
    plt.title("Mean energy for temperature=1 [kT/J] (Disordered initial state)")

    plt.figure()
    plt.plot(mc_steps, mag_temp1_random)
    plt.title("Mean magnetisation for temperature=1 [kT/J] (Disordered initial state)")

    plt.figure()
    plt.plot(mc_steps, energy_temp24_random)
    plt.title("Mean energy for temperature=2.4 [kT/J] (Disordered initial state)")

    plt.figure()
    plt.plot(mc_steps, mag_temp24_random)
    plt.title("Mean magnetisation for temperature=2.4 [kT/J] (Disordered initial state)")
    
    plt.figure()
    plt.plot(mc_steps, energy_temp1_not_random)
    plt.title("Mean energy for temperature=1 [kT/J] (Ordered initial state)")

    plt.figure()
    plt.plot(mc_steps, mag_temp1_not_random)
    plt.title("Mean magnetisation for temperature=1 [kT/J] (Ordered initial state)")

    plt.figure()
    plt.plot(mc_steps, energy_temp24_not_random)
    plt.title("Mean energy for temperature=2.4 [kT/J] (Ordered initial state)")

    plt.figure()
    plt.plot(mc_steps, mag_temp24_not_random)
    plt.title("Mean magnetisation for temperature=2.4 [kT/J] (Ordered initial state)")
    plt.show()



#calling
#plot_accepted_temp("accepted_temp_random.txt", "accepted_temp_not_random.txt")

#plot_accepted_mc("accepted_mc_random_temp1.txt", \
#                 "accepted_mc_random_temp24.txt", \
#                 "accepted_mc_not_random_temp1.txt", \
#                 "accepted_mc_not_random_temp24.txt")

plot_mean("4C_temp1_random.txt", "4C_temp24_random.txt", \
                 "4C_temp1_not_random.txt", "4C_temp24_not_random.txt")
