import numpy as np
import matplotlib.pyplot as plt

# User-defined parameters
kappa = 0.41
B = 5.2
u_tau = 6.37309e-02
Re_tau = 182.088

# Extract the data
y, U, k, eps, v2, f = np.genfromtxt("../data/test/init.dat",
                                    usecols=[0,1,2,3,4,5], unpack=True)
yp = np.multiply(y,Re_tau)

# Build the plot
end = 0.3/Re_tau  # end of the viscous low law region
U_powerlaw = yp
U_loglaw = np.divide(np.log(yp),kappa)+B
plt.semilogx(yp, U, 'k', label="Data")
plt.semilogx(yp, U_powerlaw, 'r:', label="$u^+ = y^+$")
plt.semilogx(yp, U_loglaw, 'b:', label="$u^+ = 1/\kappa \ln y^+ + B$")
plt.axis([np.min(yp),np.max(yp),0,25])

# Plot k, eps, v2, f
fig, axarr = plt.subplots(2,2)
axarr[0,0].semilogx(yp, U, label='$U$')
axarr[0,1].plot(yp, k, label='$k$')
axarr[0,1].plot(yp, v2, label='$\overline{v^2}$')
axarr[0,0.25].plot(yp, eps, label='$\epsilon$')
axarr[1,1].plot(yp, np.multiply(k,f), label='$kf$')
plt.show()

plt.show()
