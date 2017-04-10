import numpy as np
import matplotlib.pyplot as plt

# User-defined parameters
k = 0.41
B = 5.2
u_tau = 4.58794e-02
Re_tau = 1994.756

# Extract the data
y, U = np.genfromtxt("../data/init.dat", usecols=[0,1], unpack=True)
yp = np.multiply(y,Re_tau)

# Build the plot
end = 0.3/Re_tau  # end of the viscous low law region
U_powerlaw = yp
U_loglaw = np.divide(np.log(yp),k)+B
plt.semilogx(yp, U, 'k', label="Data")
plt.semilogx(yp, U_powerlaw, 'r:', label="$u^+ = y^+$")
plt.semilogx(yp, U_loglaw, 'b:', label="$u^+ = 1/\kappa \ln y^+ + B$")
plt.axis([np.min(yp),np.max(yp),0,25])

plt.show()
