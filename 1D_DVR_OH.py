import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

H_mass = 1.007825
O_mass = 15.994915
h_bar = 1.0
ang2bohr = (1.e-10)/(5.291772106712e-11)
Avo_num = 6.0221367e23
me = 9.10938356e-31
har2wave = 219474.6
c = 2.99792458e10
m_red = ((H_mass*O_mass)/(H_mass+O_mass))/(Avo_num*1000)
k = ((3600.*2.*np.pi*c)**2)*(m_red)
print "goal = " + str(np.sqrt(k/(m_red))/(2.*np.pi*c))


def grid(a, b, N):
    a = a*ang2bohr
    b = b*ang2bohr
    x = np.linspace(a,b,num=N)
    return x


def Potential_HO(f, grid):
    x = grid
    k = (f/me)*(2.41884326505e-17**2)  # conversion from seconds to atomic units
    pot = (1./2.)*k*(x**2)
    v_final = np.diag(pot)
    return v_final


def Kinetic_Calc(m1, m2, grid):
    mass_1 = m1/(Avo_num*me*1000)
    mass_2 = m2/(Avo_num*me*1000)
    m_red = (mass_1*mass_2)/(mass_1+mass_2)
    a = grid[0]
    b = grid[-1]
    N = len(grid)
    coeff = (h_bar**2/((2.*m_red)/(((float(N)-1.)/(b-a))**2)))

    Tii = np.zeros(N)

    Tii += coeff*((np.pi**2.)/3.)
    T_initial = np.diag(Tii)
    for i in range(1,N):
        for j in range(i):
            T_initial[i, j] = coeff*((-1.)**(i-j))*(2./((i-j)**2))
    T_final = T_initial + T_initial.T - np.diag(Tii)
    return T_final


def Energy(T, V):
    H = (T + V)
    En, Eigv = np.linalg.eigh(H)
    print "Ground State = " + str((En[0]) * har2wave)
    print "result = " + str((En[1]-En[0]) * har2wave)
    return En, Eigv


def Plot_wavefn(En, Eigv, V, grid):
    Eigv += En
    Eigv = Eigv*har2wave  # for scale
    Eo = np.zeros(1000)
    Eo += (En[0])*har2wave
    y = np.diag(V)*har2wave
    x = grid/ang2bohr
    plt.plot(x, y, "r-", x, Eigv[:,0], "b-", x, Eo, "k-")
    plt.xlabel("R (Angstroms)")
    plt.ylabel("Energy (cm$^{-1}$)")
    plt.show()
    plt.close()


def run():
    g = grid(-1, 1, 1000)
    V = Potential_HO(k, g)
    T = Kinetic_Calc(H_mass, O_mass, g)
    En, Eig = Energy(T, V)
    Plot_wavefn(En, Eig, V, g)


run()

# Energy()
# Plot_wavefn()
# Potential_HO(725,-1.5,1.5,100)
