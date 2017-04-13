import numpy as np
import matplotlib.pyplot as plt

# User-defined parameters
kappa = 0.384
B = 4.27

Re_tau = 5200
if (Re_tau == 5200):
    u_tau = 4.14872e-02
elif (Re_tau == 180):
    u_tau = 6.37309e-02

class vars:
    def __init__(self, filename):
        self.y, self.U, self.k, self.eps, self.v2, self.f = \
                np.genfromtxt(filename, usecols=[0,1,2,3,4,5], unpack=True)
        self.eps /= Re_tau
        self.yp = np.multiply(self.y, Re_tau)

    def plot_BL(self, fig, ax):
        end = 0.35*Re_tau  # end of the viscous low law region
        U_powerlaw = self.yp
        U_loglaw = np.divide(np.log(self.yp),kappa)+B
        ax.semilogx(self.yp, self.U, 'k', label="$\overline{v^2}-f$ results")
        ax.semilogx(self.yp, U_powerlaw, 'r:', label="$u^+ = y^+$")
        ax.semilogx(self.yp, U_loglaw, 'b:', label="$u^+ = \\frac{1}{\kappa} \ln y^+ + B$")
        ax.set_ylabel('$u^+$')
        ax.set_xlabel('$y^+$')
        ax.axis([np.min(self.yp),end,0,30])

    def plot_U(self, fig, ax, symbol, postfix):
        ax.semilogx(self.yp, self.U, symbol, label='$u^+$, '+postfix)
        ax.set_ylabel('$u^+$')
        ax.set_xlabel('$y^+$')        
                
    def plot_k_y(self, fig, ax, symbol, postfix):
        ax.plot(self.y, self.k, symbol+'--', label='$k$, '+postfix)
        ax.plot(self.y, self.v2, symbol+'-', label='$\overline{v^2}$, '+postfix)
        ax.set_ylabel('$k$, $\overline{v^2}$')
        ax.set_xlabel('$\eta$')

    def plot_k_yp(self, fig, ax, symbol, postfix):
        ax.semilogx(self.yp, self.k, symbol+'--', label='$k$, '+postfix)
        ax.semilogx(self.yp, self.v2, symbol+'-', label='$\overline{v^2}$, '+postfix)
        ax.set_ylabel('$k$, $\overline{v^2}$')
        ax.set_xlabel('$y^+$')

    def plot_eps(self, fig, ax, symbol, postfix):
        ax.semilogx(self.yp, self.eps, symbol, label='$\\varepsilon^+$, '+postfix)
        ax.set_ylabel('$\\varepsilon^+$')
        ax.set_xlabel('$y^+$')

# Extract the data
v2f = vars("../data/v2fResults_5200.dat")
DNS = vars("../data/init_5200.dat")

# Build the plot
fig, ax = plt.subplots()
v2f.plot_BL(fig, ax)
ax.legend(loc='best')
fig.set_size_inches(7,5, forward=True)
fig.savefig('../doc/turbulence/figs/loglaw_5200.png')
plt.show()

# Plot k, eps, v2, f
fig1, ax = plt.subplots()
v2f.plot_U(fig1, ax, 'r', '$\overline{v^2}-f$')
DNS.plot_U(fig1, ax, 'b', 'DNS')
ax.legend(loc='best')
fig1.tight_layout()
fig1.set_size_inches(7,5, forward=True)
fig1.savefig('../doc/turbulence/figs/U_5200.png')
plt.show()

fig2, ax = plt.subplots()
v2f.plot_k_y(fig2, ax, 'r', '$\overline{v^2}-f$')
DNS.plot_k_y(fig2, ax, 'b', 'DNS')
ax.legend(loc='best')
fig2.tight_layout()
fig2.set_size_inches(7,5, forward=True)
fig2.savefig('../doc/turbulence/figs/k_y_5200.png')
plt.show()

fig3, ax = plt.subplots()
v2f.plot_k_yp(fig2, ax, 'r', '$\overline{v^2}-f$')
DNS.plot_k_yp(fig2, ax, 'b', 'DNS')
ax.legend(loc='best')
fig3.tight_layout()
fig3.set_size_inches(7,5, forward=True)
fig3.savefig('../doc/turbulence/figs/k_yp_5200.png')
plt.show()

fig4, ax = plt.subplots()
v2f.plot_eps(fig2, ax, 'r', '$\overline{v^2}-f$')
DNS.plot_eps(fig2, ax, 'b', 'DNS')
ax.legend(loc='best')
ax.set_ylim([0,np.max(DNS.eps)*1.1])
fig4.tight_layout()
fig4.set_size_inches(7,5, forward=True)
fig4.savefig('../doc/turbulence/figs/eps_5200.png')
plt.show()
