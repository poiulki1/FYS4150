import telnetlib
import numpy as np
import sys



class planetInfo:

    def __init__(self):
        pass


    def printInfoToFile(self,startdate,enddate,IDs=[],filename = ""):



        if len(IDs) == 0 and len(filename) >0:
            IDs = np.loadtxt(filename)
        elif len(IDs) == 0 and len(filename) == 0:
            print "No IDs or filename given"
            sys.exit(1)
        print "Getting info for %g planet(s)" %len(IDs)


        kg_to_MS = 5.02739933e-31

        data = np.zeros((len(IDs),2,3))
        names = []
        masses = []
        IDs = np.array(IDs)
        for i in range(IDs.size):
            name,mass,data[i,:,:] = self.getInfo(int(IDs[i]),startdate,enddate)
            names.append(name)
            masses.append(mass)
            print "Found info for %s" %name


        #Normalizes all the masses with the current mass of the sun(if the suns mass is retrived from NASA)
        if (names[0] == "Sun"):
            kg_to_MS = 1.0/(masses[0])
        data_ = data.reshape(len(names),6)

        with open("planetData.txt","w") as outFile:


            for i in range(len(names)):

                outFile.write(str(names[i]) + "  ")
                for j in range(6):
                    outFile.write(str(data_[i,j]) + " ")
                outFile.write(str(masses[i]*kg_to_MS) + "\n")


        with open("planetConfig_" + str(len(names)) + ".txt","w") as outFile:

            outFile.write(str(len(names)) + "\n")
            for name in names:
                outFile.write(name + "\n")



    def getInfo(self,ID,startdate,enddate):

        self.t = telnetlib.Telnet()
        self.t.open('horizons.jpl.nasa.gov', 6775)

        expect = ( ( r'Horizons>', str(ID) + '\n' ),
           ( r'< Scroll.*:', 'q' ),
           ( r'Select.*E.phemeris.*:', 'E\n'),
           ( r'Observe.*:', 'v\n' ),
           ( r'Coordinate center.*:', '500@0\n' ),
           ( r'Confirm selected station.*>', 'y\n'),
           ( r'Reference plane.*:', 'eclip\n'),
           ( r'Starting .* :', str(startdate)+'\n' ),
           ( r'Ending .* :', str(enddate)+'\n' ),
           ( r'Output interval.*:', '1d\n' ),
           ( r'Accept default output.* :', 'y\n' ),
           ( r'Scroll . Page: .*%', ' '),
           ( r'Select\.\.\. .A.gain.* :', 'x\n' ))




        #Everything from the telnet site is dumped to a file
        with open('temp.txt', 'w') as fp:
            while True:
                try:
                    answer = self.t.expect(list(i[0] for i in expect), 10)
                except EOFError:
                    break


                fp.write(answer[2])
                fp.flush()
                self.t.write(expect[answer[0]][1])

        data = np.zeros((2,3))
        name = ""
        exponant = 0
        mass = 0
        foundMass = False

        #This loop looks through the dump file for the relevant data
        with open('temp.txt', 'r') as fp:
            for line in fp:
                words = line.split()
                if len(words) == 0:
                    continue


                #Look for the name of the body
                if words[0] == "Target":
                    name = words[3]
                    continue

                #This loop looks for the mass. Due to inconsistent formating from NASA, this loop is fairly complicated
                if (("Mass" in words) or ("Mass," in words)) and (not foundMass):
                    for i in range(len(words)):
                        if (words[i] == "Mass") or (words[i] == "Mass,"):

                            try:
                                if len(words[i+1]) == 6:
                                    exponant = float(words[i+1][4:])
                                else:

                                    exponant = float(words[i+1][3:])

                                if (words[i+2] == "g") or (words[i+2] == "g)"):
                                    exponant -= 3
                            except:
                                if len(words[i+2]) == 6:
                                    exponant = float(words[i+2][4:])
                                else:

                                    exponant = float(words[i+2][3:])

                                if (words[i+2] == "g") or (words[i+2] == "g)"):
                                    exponant -= 3
                            for j in range(len(words)-i):
                                try:
                                    mass = float(words[j+i])*10**exponant
                                    break
                                except:
                                    try:
                                        nwords = words[j+i].split("+")

                                        mass = float(nwords[0])*10**exponant
                                        break
                                    except:
                                        continue
                            foundMass = True
                            break


                #The initial position and velocity is always marked with "$$SOE"
                if words[0] == "$$SOE":
                    fp.next()
                    data[0,:] = np.array(fp.next().split())
                    data[1,:] = np.array(fp.next().split())
                    data[1,:] *= 365.24
                    break

        return name,mass,data


info = planetInfo()


#The info is always retrieved from the startdate, the end date is only given for the telnet site to work
start = "2018-Oct-24"
end = "2018-Nov-23"


IDs = [10,199, 299, 399, 499, 599, 699, 799,899, 999]
info.printInfoToFile(start,end, IDs, filename = "planetinfo.txt")
