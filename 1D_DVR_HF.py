import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def Potential_HO(k, a, b, N):
    x = np.linspace(a,b,num=N)
    k = k/(9.10938970000e-31)*(2.41884326505e-17)**2
    pot = (1./2.)*k*(x**2)
    v_final = np.diag(pot)
    return v_final, pot, x

def Kinetic_Calc(N, m1, m2, a, b):
    h_bar = 1.0
    mass_1 = m1/(6.02213670000e23*9.10938970000e-28)
    mass_2 = m2/(6.02213670000e23*9.10938970000e-28)
    m_red = (mass_1*mass_2)/(mass_1+mass_2)
    coeff = (h_bar)**2.0/(2.0*m_red)*(1.0/(b-a)**2.0)*(N**2)

    Tii = np.zeros(N)

    Tii += coeff*(np.pi**2./3.)
    T_initial = np.diag(Tii)
    for i in range(1,N):
        for j in range(i):
            T_initial[i,j] = coeff*(-1.0)**(i-j)*(2./(i-j)**2)
    T_final = T_initial + T_initial.T -np.diag(Tii)
    return T_final

def Energy():
    Pot, y, x = Potential_HO(725,-1,1,100)
    Kinetic = Kinetic_Calc(100, 1.007825, 15.994915, -1, 1)
    H = (Kinetic + Pot) #* 219474.6
    En, Eigv = np.linalg.eigh(H)
    Eigv = (-1.)*Eigv
    # plt.plot(x, Eigv[:,0])
    # plt.show()
    # plt.close()
    # print En[1]-En[0]
    return Eigv, En, Pot, y, x

def Plot_wavefn():
    Eigv, En, Pot, y, x = Energy()
    Eigv += En
    # V_min = min(y)
    # y -= V_min
    y = y*219474.6
    # Eigv -= V_min
    Eigv = Eigv*219474.6
    Eo = np.zeros(100)
    Eo += (En[0])*219474.6
    E1 = np.zeros(100)
    E1 += (En[1]) * 219474.6
    plt.plot(x, y, "r-", x, Eigv[:,1], "b-", x, E1, "k-")
    # plt.legend(["Potential", "Ground State Wavefunciton", "Ground State Energy"])
    plt.xlabel("R (Angstoms)")
    plt.ylabel("Energy (cm$^{-1}$)")
    plt.show()
    plt.close()

# Energy()
Plot_wavefn()