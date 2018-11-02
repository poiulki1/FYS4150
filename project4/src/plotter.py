import numpy as np
import matplotlib.pyplot as plt

def plot_mean(filename1, filename2):
    mean_temp1 = np.loadtxt(filename1)
    mean_temp24 = np.loadtxt(filename2)

    mc_steps = mean_temp1[:,0]

    energy_temp1 = mean_temp1[:,1]
    mag_temp1 = mean_temp1[:,2]

    energy_temp24 = mean_temp24[:,1]
    mag_temp24 = mean_temp24[:,2]

    plt.figure()
    plt.plot(mc_steps, energy_temp1)
    plt.title("Mean energy for temperature=1 [kT/J]")

    plt.figure()
    plt.plot(mc_steps, mag_temp1)
    plt.title("Mean magnetisation for temperature=1 [kT/J]")

    plt.figure()
    plt.plot(mc_steps, energy_temp24)
    plt.title("Mean energy for temperature=2.4 [kT/J]")

    plt.figure()
    plt.plot(mc_steps, mag_temp24)
    plt.title("Mean magnetisation for temperature=1 [kT/J]")
    plt.show()
    
plot_mean("4C_temp1_random.txt", "4C_temp24_random.txt")
