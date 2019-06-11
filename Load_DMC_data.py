import numpy as np
import matplotlib.pyplot as plt

def load_Psi(file_name_Psi, file_name_Energy, file_name_pot):
    Psi = np.load("%s.npy" %file_name_Psi)
    Energy = np.load("%s.npy" % file_name_Energy)
    Pot = np.load("%s.npy" % file_name_pot)
    return Psi, Energy, Pot


Psi, Energy, Pot = load_Psi("DMC_HO_water_Psi", "DMC_HO_water_Energy", "DMC_HO_water_pot")

me = 9.10938356e-31
Avo_num = 6.0221367e23
m_O = 15.994915 / (Avo_num*me*1000)
m_H = 1.007825 / (Avo_num*me*1000)
m_red = (m_O*m_H)/(m_O+m_H)
c = 2.99792458e10
k = (((3600.*2.*np.pi*c)**2)*m_red)*2.41884326505e-17**2
# alpha = np.sqrt(m_red*k)
har2wave = 219474.6
omega = 3600./har2wave
alpha = m_red*omega
blah = (alpha/np.pi)**(1./4.)

class wvfn(object):
    Psi = np.zeros(4)

    def __init__(self, Psi):
        self.walkers = Psi[0, :]
        self.x = Psi[1, :]
        self.weights = Psi[2, :]
        self.d = Psi[3, :]


def make_psi(Psi):
    psi = wvfn(Psi)
    return psi


psi = make_psi(Psi)

xsrqd = sum(psi.d * psi.x**2)/sum(psi.d)

N = 1./len(psi.x)
x = np.arange(-1., 1.+N, N)
analytic = np.exp(-alpha*x**2./2.)*(alpha/np.pi)**(1./4.)
Psi_sqrd = analytic**2.
ampansqrd, xxansqrd = np.histogram(x, weights=Psi_sqrd, bins=75, range=(-1., 1.), density=True)
ampan, xxan = np.histogram(x, weights=analytic, bins=75, range=(-1., 1.), density=True)
amp, xxp = np.histogram(psi.x, weights=psi.weights, bins=75, range=(-1., 1.), density=True)
ampsqrd, xxpsqrd= np.histogram(psi.x, weights=psi.d, bins=75, range=(-1., 1.), density=True)
binsp = (xxp[1:]+xxp[:-1])/2.
binspsqrd = (xxpsqrd[1:]+xxpsqrd[:-1])/2.
binsan = (xxan[1:]+xxan[:-1])/2.
binsansqrd = (xxansqrd[1:]+xxansqrd[:-1])/2.

plt.plot(binsan, ampan, 'k', label='Psi')
plt.plot(binsansqrd, ampansqrd, 'b', label='Psi squared')
plt.plot(binsp, amp, 'g', label='Psi DMC')
plt.plot(binspsqrd, ampsqrd, 'r', label='Psi DMC squared')
plt.legend()
plt.show()
print np.mean(Energy[500:]*har2wave)
print (1./2.)*alpha**(-1.)
print xsrqd

