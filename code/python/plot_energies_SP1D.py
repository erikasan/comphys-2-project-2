import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

omegas = np.linspace(1, 20, 5)
hidden_layers = [1]#[1, 10, 100]

sns.set()

exact = omegas/2
plt.plot(omegas, exact, label="Analytical solution")

for h in hidden_layers:
    energy = np.loadtxt("../../output/energies_SP1D_imp_h={}".format(h))
    plt.plot(omegas, energy, label="{} hidden layers".format(h))

plt.xlabel(r'$\omega$')
plt.ylabel(r'$E_0$')
plt.legend()
plt.show()
