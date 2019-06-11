import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

H_mass = 1.007825
O_mass = 15.994915
h_bar = 1
ang2bohr = 1.88973
Ava_num = 6.02213670000e23
me = 9.10938970000e-31
har2wave = 219474.6

def Potential_Calc(file_name, N, a, b):
    x = []
    y = []
    with open("%s" %file_name, "r") as f:
        for line in f:
            words = line.split()
            x.append(float(words[0]))
            y.append(float(words[1]))

    int = interpolate.splrep(x,y,s=0)
    xnew = np.linspace(a, b, num=N)
    ynew = interpolate.splev(xnew, int, der=0)
    index = np.argmin(ynew)

    print xnew[index]
    v_final = np.diag(ynew)
    return v_final, ynew, xnew

def Kinetic_Calc(steps, a, b):
    m_O = 15.994915/(6.02213670000e23*9.10938970000e-28)
    m_H = 1.007825/(6.02213670000e23*9.10938970000e-28)
    m_OH = m_O+m_H
    m_red = (m_OH*m_H)/(m_OH+m_H)
    a = a*1.88973
    b = b*1.88973
    coeff = (h_bar)**2.0/(2.0*m_red)*(1.0/(b-a)**2.0)*(steps**2)

    Tii = np.zeros(steps)

    Tii += coeff*(np.pi**2/3)
    T_initial = np.diag(Tii)
    for i in range(1,steps):
        for j in range(i):
            T_initial[i,j] = coeff*(-1.0)**(i-j)*(2./(i-j)**2)
    T_final = T_initial + T_initial.T -np.diag(Tii)
    return T_final

def Energy():
    Pot, ynew, xnew = Potential_Calc("SPE_W_data.txt", 1000, 0.5, 1.5)
    Kinetic = Kinetic_Calc(1000, 0.5, 1.5)
    H = (Kinetic + Pot)
    En, Eigv = np.linalg.eigh(H)
    print En[0]
    print (En[1]-En[0])*har2wave
    return Eigv, En, ynew, xnew

def Plot_wavefn():
    Eigenv, En, y, x = Energy()
    Eigenv += En
    V_min = min(y)
    y -= V_min
    y = y*219474.6
    Eigenv -= V_min
    Eigenv = Eigenv*219474.6
    Eo = np.zeros(1000)
    Eo += (En[0] -V_min)*219474.6
    plt.plot(x, y, "r-", x, Eigenv[:,0], "b-", x, Eo, "k-")
    plt.legend(["Potential", "Ground State Wavefunciton", "Ground State Energy"])
    plt.xlabel("R (Angstoms)")
    plt.ylabel("Energy (cm$^{-1}$)")
    plt.show()
    plt.close()

Plot_wavefn()


# Energy()
# Kinetic_Calc(500, 0.7, 1.3), Potential_Calc("SPE_W_data.txt", 500, 0.7, 1.3)
# Potential_Calc("SPE_W_data.txt")
# Kinetic_Calc(100,0.7,1.3)
