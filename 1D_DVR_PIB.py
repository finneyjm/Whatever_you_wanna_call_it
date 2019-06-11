import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def Potential_PIB(l, a, b, N):
    x = np.linspace(a,b,num=N)
    l = np.linspace(-l/2., l/2., num=l*N/((b-a)/2))

    v_final = np.diag(pot)

    # plt.plot(x, pot, "b-")
    # plt.show()
    # plt.close()
    return v_final, pot, x

def Kinetic_Calc(m1, m2, a, b, N):
    h_bar = 1.0
    mass_1 = m1/(6.02213670000e23*9.10938970000e-28)
    mass_2 = m2/(6.02213670000e23*9.10938970000e-28)
    m_red = (mass_1*mass_2)/(mass_1+mass_2)
    a = a*1.88973
    b = b*1.88973
    print b-a
    coeff = (h_bar)**2.0/(2.0*m_red)*(1.0/(b-a)**2.0)*(float(N)**2)

    Tii = np.zeros(N)

    Tii += coeff*((np.pi**2.)/3.)
    T_initial = np.diag(Tii)
    for i in range(1,N):
        for j in range(i):
            T_initial[i,j] = coeff*((-1.0)**(i-j))*(2./((i-j)**2))
    T_final = T_initial + T_initial.T -np.diag(Tii)
    return T_final

def Energy():
    # Pot, y, x = Potential_HO(725., -1.5, 1.5, 1000)
    Kinetic = Kinetic_Calc(1.007825, 15.994915, -1.5, 1.5, 100)
    H = (Kinetic+ Pot) #* 219474.6
    En, Eigv = np.linalg.eigh(H)
    # Eigv = (-1.)*Eigv
    # plt.plot(x, Eigv[:,0])
    # plt.show()
    # plt.close()
    print En[1]
    return Eigv, En, #Pot, y, x

def Plot_wavefn():
    Eigv, En= Energy()
    Eigv += En
    # V_min = min(y)
    # y -= V_min
    y = y*219474.6
    # Eigv -= V_min
    Eo = np.zeros(1000)
    Eo += (En[0])*219474.6
    # E1 = np.zeros(100)
    # E1 += (En[1]) * 219474.6
    plt.plot(x, y, "r-", x, Eigv[:,0], "b-", x, Eo, "k-")
    # plt.legend(["Potential", "Ground State Wavefunciton", "Ground State Energy"])
    plt.xlabel("R (Angstoms)")
    plt.ylabel("Energy (cm$^{-1}$)")
    plt.show()
    plt.close()

Energy()