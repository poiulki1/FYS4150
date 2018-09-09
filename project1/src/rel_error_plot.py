import numpy as np
import matplotlib.pyplot as plt
import seaborn

def read_2_col_file(file_name):
    """
    This function reads the text file line wise and splits the two columns and
    stores it in separate lists. We want to plot 1/h, so the column two is
    implemented to calculate this already in the read stage.
    """
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
    """
    Logarithmic plot of the data read from the text file stored in two lists,
    file name needs to be given here in the code.
    """
    error_LU, h_LU = read_2_col_file("error_LU_data.txt") 
    error, h = read_2_col_file("error_data.txt")
    
    log_h_LU = np.log10(h_LU)
    log_error_LU = np.log10(error_LU)
    log_h = np.log10(h)
    log_error = np.log10(error)

    plt.plot(log_h_LU, log_error_LU, label="LU-decomposition")
    plt.plot(log_h, log_error, label="general Thomas algo.")
    plt.legend()
    plt.title("Logarithmic plot of relative error ")
    plt.xlabel("$log_{10}(time \, step)$")
    plt.ylabel("$log_{10}(error)$")
    plt.show()


main()

