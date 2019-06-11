import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import copy
from scipy.optimize import curve_fit

# DMC parameters
dtau = 5.
N_0 = 5000
time_total = 10000.
alpha = 1./(2.*dtau)

# parameters for the potential
me = 9.10938356e-31
Avo_num = 6.0221367e23
m_O = 15.994915 / (Avo_num*me*1000)
m_H = 1.007825 / (Avo_num*me*1000)
m_red = (m_O*m_H)/(m_O+m_H)
c = 2.99792458e10
k = (((3600.*2.*np.pi*c)**2)*m_red)*2.41884326505e-17**2

# factor for random walk
sigma = np.sqrt(dtau/m_red)

# conversion factor
har2wave = 219474.6


# Make the object that contains the starting wave function
class Psi(object):
    walkers = 0

    def __init__(self, walkers):
        self.walkers = np.linspace(0, walkers-1, num=walkers)
        self.x = np.zeros(walkers) + 0.96
        self.weights = np.zeros(walkers) + 1.
        self.d = np.zeros(walkers)
        self.weights_i = np.zeros(walkers) + 1.


def make_Psi(walkers):
    psi = Psi(walkers)
    return psi


wvfn = make_Psi(N_0)


# Kinetic energy is a random walk of the walker along a Gaussian in the x direction
def Kinetic(Psi):
    N_t = len(Psi.x)
    random_walk = np.random.normal(0.0, sigma, N_t)
    Psi.x += random_walk
    return Psi


# Pulls the potential from a given file built from my Gaussian_output function
def Build_Potential(file_name):
    x = []
    y = []
    with open("%s" %file_name, "r") as f:
        for line in f:
            words = line.split()
            x.append(float(words[0]))
            y.append(float(words[1]))
    return x, y


def morse(x, De, a, re, v):
    return De*(np.exp(-2.*a*(x-re)) - 2*np.exp(-a*(x-re))) + v


# Calculates the potential for each walker along the interpolated potential from the points built from the last function
def Potential_calc(x, y, Psi):
    Psi.x = np.array(Psi.x, dtype='float64')  # converts the x values into the right type for the input of the

    paras, covar = curve_fit(morse, x, y)

    Vi = morse(Psi.x, paras[0], paras[1], paras[2], paras[3])

    # int = interpolate.splrep(x, y, s=0)
    # Vi = interpolate.splev(Psi.x, int, der=0)
    return Vi


# Calculate Vref/Er using Er = mean(V) -alpha*(P-P0)/P0
def V_ref_calc(Vi,Psi):
    P0 = sum(Psi.weights_i)
    P = sum(Psi.weights)
    V_ref = sum(Psi.weights*Vi)/P - alpha*(sum((Psi.weights-Psi.weights_i))/P0)
    return V_ref


# Calculating the weights of the walkers
def Potential(Vi, Vref, Psi):
    Psi.weights = Psi.weights*np.exp(-(Vi-Vref)*dtau)
    # Condition to prevent one walker from obtaining all the weight
    threshold = 1/sum(Psi.weights_i)
    death = np.argwhere(Psi.weights < threshold)
    for i in death:
        ind = np.argmax(Psi.weights)
        Biggo_weight = copy.deepcopy(Psi.weights[ind])
        Biggo_pos = copy.deepcopy(Psi.x[ind])
        Psi.weights[i[0]] = Biggo_weight/2.
        Psi.weights[ind] = Biggo_weight/2.
        Psi.x[i[0]] = Biggo_pos
    return Psi


# A similar weight calculation that also tracks the descendants of the walkers that are being replicated
def dPot(Vi, Vref, Psi):
    Psi.weights = Psi.weights*np.exp(-(Vi-Vref)*dtau)
    # Condition to prevent one walker from obtaining all the weight
    threshold = 1/sum(Psi.weights_i)
    death = np.argwhere(Psi.weights < threshold)
    for i in death:
        ind = np.argmax(Psi.weights)
        Biggo_weight = copy.deepcopy(Psi.weights[ind])
        Biggo_pos = copy.deepcopy(Psi.x[ind])
        Biggo_num = copy.deepcopy(Psi.walkers[ind])
        Psi.weights[i[0]] = Biggo_weight/2.
        Psi.weights[ind] = Biggo_weight/2.
        Psi.walkers[i] = Biggo_num
        Psi.x[i[0]] = Biggo_pos
    return Psi


# Creates a new value for the Psi object that give descendant weights to the walkers after some dtau
def descendants(Psi):
    for i in xrange(len(Psi.walkers)):
        Psi.d[i] = np.sum(Psi.weights[Psi.walkers == i])
    return Psi.d


# The function that runs everything
def run(propagation):
    Psi = Kinetic(wvfn)
    x, y = Build_Potential("SPE_W_data.txt")  # pulls the values from the txt file to be continuously used
    Vi = Potential_calc(x, y, Psi)
    Eref = np.array([])
    Vref = V_ref_calc(Vi, Psi)
    Eref = np.append(Eref, Vref)
    new_psi = Potential(Vi, Vref, Psi)

    # initial parameters before running the calculation
    DW = False  # a parameter that when True will implement descendant weighting
    Psi_dtau = 0
    for i in xrange(int(time_total)):
        if DW is False:
            prop = float(propagation)

        Psi = Kinetic(new_psi)
        Vref = V_ref_calc(Vi, Psi)
        Vi = Potential_calc(x, y, Psi)
        Eref = np.append(Eref, Vref)

        if DW is False:
            new_psi = Potential(Vi, Vref, Psi)
        elif DW is True:
            if Psi_dtau == 0:
                Psi_tau = copy.deepcopy(Psi)
                Psi_dtau = copy.deepcopy(Psi_tau)
                new_psi = dPot(Vi, Vref, Psi_dtau)
            else:
                new_psi = dPot(Vi, Vref, Psi)
            prop -= dtau

        if i >= (time_total - float(propagation)) and prop > 0:  # start of descendant weighting
            DW = True
        elif i >= (time_total - float(propagation)) and prop == 0:  # end of descendant weighting
            d_values = descendants(new_psi)
            Psi_tau.d += d_values

    psi = np.zeros((4, N_0))
    psi[0:] = Psi_tau.walkers
    psi[1:] = Psi_tau.x
    psi[2:] = Psi_tau.weights
    psi[3:] = Psi_tau.d
    np.save("DMC_1dpot_water_Psi", psi)
    np.save("DMC_1dpot_water_Pot", Vi)
    np.save("DMC_1dpot_water_Energy", Eref)
    return


run(50)

