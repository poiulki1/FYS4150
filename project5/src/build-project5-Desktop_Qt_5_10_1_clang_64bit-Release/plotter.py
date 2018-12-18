import matplotlib.pyplot as plt
import numpy as np

print("do you wish to plot rk4 and mc or variance? ")
argument = str(input("mc and rk4: yes, if no variance is plotted "))

if argument == "yes":

    #arg1 = str(input("Is this monte carlo simulation or runge kutta? type mc or rk4 "))
    number_arg = str(input("What is the value used for a,c,d,di,e,f ? "))
    number_arg1 = number_arg.split("_")


    something1 = []
    something2 = []
    something3 = []

    for i in range(1,5):
        something1.append(np.genfromtxt("../result/run_mc_b_" + str(i) + "_acddief_"+ number_arg +  ".txt"))
        something2.append(np.genfromtxt("../result/run_rk4_b_" + str(i) + "_acddief_"+ number_arg +  ".txt"))
        #something3.append(np.genfromtxt("../result/sigma_mc_b_" + str(i) + "_acddief_"+ number_arg +  ".txt"))

    for i in range(1,5):
        plt.clf()
        plt.subplot(2,1,1)
        plt.plot(something1[i-1][:,0], something1[i-1][:,1])
        plt.plot(something1[i-1][:,0], something1[i-1][:,2])
        plt.plot(something1[i-1][:,0], something1[i-1][:,3])
        plt.legend(['Susceptible', 'infected', 'recovered'], loc = 'upper right')

        plt.title("Monte carlo: $a =$ " + str(number_arg1[0]) + " ,$b =$ " + str(i)
         + ",$ c =$ " + str(number_arg1[1])+  " , $d =$ " + str(number_arg1[2]) + " ,$d_{i} = $" + str(number_arg1[3]) +
        ", $e = $" +  str(number_arg1[4]) + " ,$f = $" + str(number_arg1[5]))

        plt.grid('on')
        plt.xlabel("time")
        plt.ylabel("Number of people")
        plt.tight_layout()


        plt.subplot(2,1,2)
        plt.plot(something2[i-1][:,0], something2[i-1][:,1])
        plt.plot(something2[i-1][:,0], something2[i-1][:,2])
        plt.plot(something2[i-1][:,0], something2[i-1][:,3])
        plt.legend(['Susceptible', 'infected', 'recovered'], loc = 'upper right')

        plt.title("Runge kutta 4: $a =$ " + str(number_arg1[0]) + " ,$b =$ " + str(i)
         + ",$ c =$ " +  str(number_arg1[1]) + " , $d =$ " + str(number_arg1[2]) + "$, d_{i} = $" + str(number_arg1[3]) +
        ", $e = $" +  str(number_arg1[4]) + " ,$f = $" + str(number_arg1[5]))
        plt.grid('on')
        plt.xlabel("time")
        plt.ylabel("Number of people")
        plt.tight_layout()
        plt.savefig("../figures/run_b_" + str(i) + "_abcddief_"+ number_arg + ".pdf", bbox_inches = "tight")
        plt.show()


elif argument == "no":
    number_arg1 = str(input("What is the value used for a,b,c,d,di,e,f ? "))
    something1 = []
    something2 = []

    for i in range(1,5):
        something1.append(np.genfromtxt("../result/variance_b" + str(i) + "_acddief_"+ number_arg1 +  ".txt"))
        something2.append(np.genfromtxt("../result/run_mc_b_" + str(i) + "_acddief_"+ number_arg1 +  ".txt"))

    for i in range(1,5):
        plt.clf()

        plt.errorbar(something1[i-1][:,0], something1[i-1][:,1], yerr = something2[i-1][:,1], marker ="s")
        plt.errorbar(something1[i-1][:,0], something1[i-1][:,2], yerr = something2[i-1][:,2], marker = "s")
        plt.errorbar(something1[i-1][:,0], something1[i-1][:,3], yerr = something2[i-1][:,3], marker = "s")

        plt.title("Monte Carlo 4: $a =$ " + str(number_arg1[0]) + " ,$b =$ " + str(i)
         + ",$ c =$ " +  str(number_arg1[1]) + " , $d =$ " + str(number_arg1[2]) + "$, d_{i} = $" + str(number_arg1[3]) +
        ", $e = $" +  str(number_arg1[4]) + " ,$f = $" + str(number_arg1[5]))
        plt.legend(['Susceptible', 'infected', 'recovered'], loc = 'best')
        #plt.show()
        plt.savefig("../figures/variance_b_" + str(i) + "_abcddief_"+ number_arg1 + ".pdf", bbox_inches = "tight")
