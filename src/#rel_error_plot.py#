import numpy as np
import matplotlib.pyplot as plt
import seaborn

def read_2_col_file(file_name):
    myfile = open(file_name, "r")
    col1 = []; col2 = []
    lines = myfile.readlines()
    for line in lines:
        value = line.split()
        col1.append(float(value[0]))
        col2.append(1/(float(value[1])))
    myfile.close()
    return col1, col2

def main():
    error, h = read_2_col_file("error_data.txt")
    log_h = np.log10(h)
    log_error = np.log10(error)

    plt.plot(log_h, log_error)
    plt.title("Logarithmic plot of relative error ")
    plt.xlabel("$log_{10}(time \, step)$")
    plt.ylabel("$log_{10}(error)$")
    plt.show()

main()

