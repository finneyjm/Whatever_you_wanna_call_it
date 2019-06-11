import numpy as np
import matplotlib.pyplot as plt

dtau = 5.
N_0 = 5000
time_total = 1000.
alpha = 1./(2.*dtau)
me = 9.10938356e-31
Avo_num = 6.0221367e23
m_O = 15.994915 / (Avo_num*me*1000)
m_H = 1.007825 / (Avo_num*me*1000)
m_red = (m_O*m_H)/(m_O+m_H)
c = 2.99792458e10
k = ((((3600.*2.*np.pi*c)**2)*(m_red))*((2.41884326505e-17)**2))
sigma = np.sqrt(dtau/m_red)

Psi_initial = np.zeros(N_0)

def Kinetic(Psi):
    N_t = len(Psi)
    random_walk = np.random.normal(0.0, sigma, N_t)
    mod_Psi = Psi + random_walk
    return mod_Psi


def Potential_HO(Psi):
    Vi = (1. / 2.) * k * (Psi ** 2)
    return Vi


def V_ref_calc(pot, Psi):
    N_t = len(Psi)
    V_ref = np.mean(pot) - (alpha * (N_t-N_0)/N_0)
    return V_ref


def Potential(Vi, Vref, Psi):
    weights = np.exp(-(Vi-Vref)*dtau)
    new_Psi = np.array([])
    for i in xrange(len(weights)):
        for n in xrange(int(weights[i])):
            new_Psi = np.append(new_Psi, Psi[i])
        bar = np.random.random()
        if bar < weights[i] - float(int(weights[i])):
            new_Psi = np.append(new_Psi, Psi[i])
    return new_Psi

def descendants(Psi):
    asdkfh

def run():
    Psi = Kinetic(Psi_initial)
    Vi = Potential_HO(Psi)
    Vref = V_ref_calc(Vi, Psi)
    new_psi = Potential(Vi, Vref, Psi)
    for i in xrange(int(time_total / dtau)):
        Psi = Kinetic(new_psi)
        Vi = Potential_HO(Psi)
        Vref = V_ref_calc(Vi, Psi)
        new_psi = Potential(Vi, Vref, Psi)
    amp, xx = np.histogram(new_psi, bins=25, range=(-0.75, 0.75), density=True)
    bins = (xx[1:]+xx[:-1])/2.
    plt.plot(bins, amp)
    print Vref
    
    plt.show()
    return new_psi


# print Psi_trial
run()

# Kinetic(Psi)
