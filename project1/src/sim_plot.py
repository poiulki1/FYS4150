import numpy as np
import matplotlib.pyplot as plt
import seaborn

def read_2_col_file(file_name):
    """
    This function is opening a given text file, reads line for line, splits
    the columns and appends to the lists. If there's a single element/column,
    it appends these elements in a different list (n_list).
    """
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
    """
    Main function, assign the read columns to the n, y and x. Then use n
    list to slice y and x for different mesh points to plot. File name needs to
    be given in the code below (not argument)
    """
    x_exact = np.linspace(0, 1, 1e5)
    n, y, x = read_2_col_file("sim_data.txt")
    legend_list = []
    start = 0; stop = n[0] + 2
    for i in range(1, len(n)):
        legend_list.append("n="+str(n[i-1]))

        plt.plot(x[start:stop], y[start:stop], label = legend_list[i-1])

        start = stop
        stop +=  2 + n[i]
        if i == (len(n)-1):
            plt.plot(x[start:stop], y[start:stop], label = "n="+str(n[-1]))

    #plotting the exact function to compare
    plt.plot(x_exact, exact(x_exact), "r--", label ="exact")
    plt.legend()
    plt.xlabel('$x$')
    plt.ylabel('$v(x)$')
    plt.title("Numerical solution of Poisson's equation with different number of mesh-points")
    plt.show()



main()
