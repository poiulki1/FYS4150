import numpy as np
import matplotlib.pyplot as plt
import seaborn

def read_2_col_file(file_name):
    myfile = open(file_name, "r")
    n_list = []
    col1 = []
    col2 = []
    lines = myfile.readlines()
    for line in lines:
        value = line.split()
        if len(value) <= 1:
            n_list.append(int(value[0]))
        else:
            col1.append(float(value[0]))
            col2.append(float(value[1]))
        myfile.close()
    return n_list, col1, col2

def exact(x):
    return 1.0 - (1 - np.exp(-10))*x -np.exp(-10*x)

def main():
    x_exact = np.linspace(0, 1, 1e5)
    n, y, x = read_2_col_file("sim_data.txt")
    legend_list = []
    start = 0; stop = n[0] + 2
    for i in range(1, len(n)):
        s = "n=", n[i-1]
        legend_list.append(s)
        
        plt.plot(x[start:stop], y[start:stop], label = legend_list[i-1])
        
        start = stop
        stop +=  2 + n[i]
        
    plt.plot(x_exact, exact(x_exact), "r--", label ="exact")
    plt.legend()
    plt.xlabel('$x$')
    plt.ylabel('$v(x)$')
    plt.title("Numerical solution of Poisson's equation with different number of mesh-points")
    plt.show()



main()
