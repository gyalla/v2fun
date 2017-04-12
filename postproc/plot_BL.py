import numpy as np
import matplotlib.pyplot as plt

# User-defined parameters
kappa = 0.41
B = 5.2
u_tau = 6.37309e-02
Re_tau = 182.088

class vars:
    def __init__(self, filename):
        self.y, self.U, self.k, self.eps, self.v2, self.f = \
                np.genfromtxt(filename, usecols=[0,1,2,3,4,5], unpack=True)
        self.yp = np.multiply(self.y, Re_tau)

    def plot_BL(self, fig, ax):
        end = 0.3/Re_tau  # end of the viscous low law region
        U_powerlaw = self.yp
        U_loglaw = np.divide(np.log(self.yp),kappa)+B
        ax.semilogx(self.yp, self.U, 'k', label="Data")
        ax.semilogx(self.yp, U_powerlaw, 'r:', label="$u^+ = y^+$")
        ax.semilogx(self.yp, U_loglaw, 'b:', label="$u^+ = 1/\kappa \ln y^+ + B$")
        ax.axis([np.min(self.yp),np.max(self.yp),0,25])

    def plot4(self, fig, axarr, symbol):
        axarr[0,0].semilogx(self.yp, self.U, symbol, label='$U$')
        axarr[0,1].plot(self.yp, self.k, symbol, label='$k$')
        axarr[0,1].plot(self.yp, self.v2, symbol, label='$\overline{v^2}$')
        axarr[1,0].plot(self.yp, self.eps/180, symbol, label='$\epsilon$')
        axarr[1,1].plot(self.yp, np.multiply(self.k,self.f)/180, symbol, 
                        label='$kf$')

# Extract the data
v2f = vars("../data/v2fResults_180.dat")
DNS = vars("../data/init.dat")

# Build the plot
fig1, ax = plt.subplots()
v2f.plot_BL(fig1, ax)
ax.legend(loc='best')
plt.show()

# Plot k, eps, v2, f
fig2, axarr = plt.subplots(2,2)
v2f.plot4(fig2, axarr, 'r--')
DNS.plot4(fig2, axarr, 'b:')
axarr[1,1].legend(['$\overline{v^2}-f$', 'DNS'],loc='best')
plt.show()
