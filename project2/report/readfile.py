import numpy as np
import matplotlib.pyplot as plt


data_repulsive = np.fromfile("eig_vec_rep_w05_max5.bin")
data_no_repulsive = np.fromfile("eig_vec_no_rep_w05_max5.bin")

j = 0


step = 200
list1 = []
list2 = []

points = len(data_repulsive)*len(data_repulsive)

for i in range(0,len(data_repulsive),step):
    list1.append(data_repulsive[j:i])
    list2.append(data_no_repulsive[j:i])
    if i == step + j:
        j += step
        continue

eigen_vectors1 = list1[1:]
eigen_vectors2 = list2[1:]

u_squared_repulsive = np.zeros(len(eigen_vectors1))
u_squared_no_repulsive = np.zeros(len(eigen_vectors2))

for k in range(len(eigen_vectors1)):
    u_squared_repulsive[k] = np.conj(eigen_vectors1[0][k])*eigen_vectors1[0][k]
    u_squared_no_repulsive[k] = np.conj(eigen_vectors2[0][k])*eigen_vectors2[0][k]



t = np.linspace(0,40, 199)
plt.title("Probability density for $\\omega = 0.001$")
plt.plot(t,u_squared_repulsive)
plt.plot(t,u_squared_no_repulsive)
plt.xlabel("scaled radial distance, $\\rho$")
plt.ylabel("Probability density, $|u(\\rho)|^{2}$")
plt.legend(["Interacting","Non-interacting"])
plt.grid('on')
plt.savefig('Probability_density_w5.pdf')
plt.show()
