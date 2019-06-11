import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

dtau = 5.
N_0 = 5000
time_total = 1000.
alpha = 1. / (2. * dtau)
me = 9.10938356e-31
Avo_num = 6.0221367e23
m_O = 15.994915 / (Avo_num * me * 1000)
m_H = 1.007825 / (Avo_num * me * 1000)
m_OH = m_H + m_O
m_red = (m_OH * m_H) / (m_OH + m_H)
sigma = np.sqrt(dtau / m_red)
har2wave = 219474.6

Psi_inital = np.zeros(N_0) + 0.96


def Kinetic(Psi):
    N_t = len(Psi)
    random_walk = np.random.normal(0., sigma, N_t)
    mod_Psi = Psi + random_walk
    return mod_Psi


def Build_Potential(file_name):
    x = []
    y = []
    with open("%s" %file_name, "r") as f:
        for line in f:
            words = line.split()
            x.append(float(words[0]))
            y.append(float(words[1]))
    return x, y


def Potential_calc(x, y, Psi):
    Psi = np.array(Psi, dtype='float64')
    int = interpolate.splrep(x, y, s=0)
    Vi = interpolate.splev(Psi, int, der=0)
    return Vi


def V_ref_calc(pot, Psi):
    N_t = len(Psi)
    V_ref = np.mean(pot) - (alpha * (N_t - N_0) / N_0)
    return V_ref


def Potential(Vi, Vref, Psi):
    weights = np.exp(-(Vi - Vref) * dtau)
    new_Psi = np.array([])
    for i in xrange(len(weights)):
        for n in xrange(int(weights[i])):
            new_Psi = np.append(new_Psi, Psi[i])
        bar = np.random.random()
        if bar < weights[i] - float(int(weights[i])):
            new_Psi = np.append(new_Psi, Psi[i])
    return new_Psi


def run():
    Psi = Kinetic(Psi_inital)
    x, y = Build_Potential("SPE_W_data.txt")
    Vi = Potential_calc(x, y, Psi)
    Eref = np.array([])
    Vref = V_ref_calc(Vi, Psi)
    Eref = np.append(Eref, Vref)
    new_psi = Potential(Vi, Vref, Psi)
    for i in xrange(int(time_total / dtau)):
        Psi = Kinetic(new_psi)
        Vref = V_ref_calc(Vi, Psi)
        Vi = Potential_calc(x, y, Psi)
        Eref = np.append(Eref, Vref)
        new_psi = Potential(Vi, Vref, Psi)
    amp, xx = np.histogram(new_psi, bins=25, range=(0.6, 1.5), density=True)
    bins = (xx[1:] + xx[:-1]) / 2.
    plt.plot(bins, amp)
    plt.show()
    plt.plot(Eref*har2wave)
    plt.show()
    E0 = np.mean(Eref[25:])
    print E0
    return new_psi


# print Psi_trial
run()
# Kinetic(Psi)
