import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

localEnergy = np.loadtxt("../../output/convergence.txt")

sns.set()
plt.plot(localEnergy)
plt.xlabel("Parameter updates")
plt.ylabel(r"$\langle E_L \rangle$")
plt.hlines(2, 0, localEnergy.size, linestyles='dashed', label=r'$E_0 = 2\omega$')
plt.legend()
plt.title("Gradient descent, two non-interacting electrons in 2D")
#plt.show()
plt.savefig('../../figures/convergence.pdf', type='pdf')
