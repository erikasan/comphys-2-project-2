import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

omegas = np.linspace(0.5, 4, 10)
hidden_layers = [1, 5, 10, 100]

sns.set()

exact = omegas/2
plt.plot(omegas, exact, label="Analytical solution")

for h in hidden_layers:
    energy = np.loadtxt("../../output/energies_SP1D_h={}".format(h))
    plt.plot(omegas, energy, label="{} hidden layers".format(h))

plt.xlabel(r'$\omega$')
plt.ylabel(r'$E_0$')
plt.legend()
plt.show()
