import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

p = 1000
t = 1.0
fname = 'test_mc/table.out'  # f'runs/{p}bar/{t:.2f}k/table.out'
data = np.loadtxt(fname, skiprows=501)
erg, alat, dens, press, xacr, vacr = data.transpose()
alat = alat/4

fig, axes = plt.subplots(3, 2, figsize=(12, 9), dpi=200, sharex=True)

axes[0, 0].plot(erg, c='#3a7')
axes[0, 0].set_ylabel('Energía por átomo (meV)')

axes[0, 1].plot(dens, c='#c33')
axes[0, 1].set_ylabel('Densidad (kg/lt)')

axes[1, 0].plot(alat, c='#37b')
axes[1, 0].set_ylabel('Parámetro de red (angs)')

axes[1, 1].plot(press, c='#649')
axes[1, 1].set_ylabel('Presión (bar)')

axes[2, 0].plot(xacr, c='#444')
axes[2, 0].set_ylabel('Tasa de aceptación desp. partícula')

axes[2, 1].plot(vacr, c='#444')
axes[2, 1].set_ylabel('Tasa de aceptación mod. volumen')

for ax in axes.flatten():
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='major', width=.8, length=4)
    ax.tick_params(which='minor', width=.6, length=3)
    ax.grid(color='grey', linestyle='-', linewidth=.25)

fig.suptitle(f'Observables NpT para p={int(p)} bar y T={t:.1f} K')
fig.supxlabel('Ciclo de MCMC')
fig.tight_layout(pad=.8)
plt.show()
