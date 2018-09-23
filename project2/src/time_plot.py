import matplotlib.pyplot as plt
import numpy as np

n = [10, 20, 50, 80, 110, 140, 170]

jacobi_time = [6.7e-5, 0.000286, 0.007528, 0.041561, 0.156877, 0.401499, 0.841207]

arma_time = [4.9e-5, 0.0001, 0.000573, 0.001519, 0.002644, 0.006798, 0.00995]

n_rotate = [10, 20, 50, 80, 110, 140, 170, 200]
jacobi_rotate = [150, 620, 4111, 10726, 20444, 33126, 48856, 67900]

plt.figure()
plt.grid("on")
plt.plot(n, jacobi_time, label="Jacobi method")
plt.plot(n, arma_time, label="Armadillo function")
plt.axis([10, 180, -0.15, 1])
plt.xlabel("Number of mesh points")
plt.ylabel("Time [s]")
plt.title("Efficiency of two algorithms to find eigenpars")
plt.legend()
#plt.savefig("timePlot.pdf")

plt.figure()
plt.grid("on")
plt.plot(n_rotate, jacobi_rotate)
plt.title("Number of Jacobi transformations as function of mesh-points")
plt.ylabel("Transformations")
plt.xlabel("Mesh-points")
#plt.savefig("transformations.pdf")
