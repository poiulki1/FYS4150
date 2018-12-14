import matplotlib.pyplot as plt
import numpy as np


arg1 = str(input("Is this monte carlo simulation or runge kutta? type mc or rk4 "))
something = []
for i in range(1,5):
    something.append(np.genfromtxt("test_" + str(i) +  ".txt"))

for i in range(1,5):
    plt.plot(something[i-1][:,0], something[i-1][:,1])
    plt.plot(something[i-1][:,0], something[i-1][:,2])
    plt.plot(something[i-1][:,0], something[i-1][:,3])
    plt.legend(['Susceptible', 'infected', 'recovered'], loc = 'best')
    plt.savefig("test_" + arg1 + "_b" + str(i) + ".pdf")
