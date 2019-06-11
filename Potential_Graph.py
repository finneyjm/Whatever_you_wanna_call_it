import matplotlib.pyplot as plt
import numpy as np

BR = []
SPE = []
BR_pot = []
Pot = []
BR_caltech = []
Pot_caltech = []

with open("SPE_W_data.txt", "r") as my_file:
    for line in my_file:
        words = line.split()
        BR.append(float(words[0]))
        SPE.append((float(words[1]) + 76.0607630898) * 219474.6)
with open("BR_pot_cut.txt", 'r') as f:
    for line1 in f:
        num = line1.split()
        BR_pot.append(float(num[0]))
        Pot.append(float(num[1]))
with open("caltech_data.txt", 'r') as f2:
    for line2 in f2:
        num2 = line2.split()
        BR_caltech.append(float(num2[0])-0.19706045305)
        Pot_caltech.append(float(num2[1]))

plt.plot(BR, SPE, 'bo', label= "MP2/aug-cc-pVTZ")
plt.plot(BR_pot, Pot, 'g-', label= "Potential")
plt.plot(BR_caltech, Pot_caltech, 'r-', label= "Caltech_modified")
plt.legend(loc = "upper right")
plt.xlabel("R (Angstroms)")
plt.ylabel("Potential Energy (cm-1)")
plt.title("Single Point Energy Calc for Water")
plt.ylim(top=30000)
plt.xlim(left=0.7)
plt.savefig("SPE_Water_new1.png")
plt.show()
plt.close()