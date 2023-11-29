import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

p = 1000
t = 1.0
fname = 'test_md/table.out'  # f'runs/{p}bar/{t:.2f}k/table.out'
data = np.loadtxt(fname, skiprows=1001)
etot, epot, ekin, temp, press = data.transpose()
time = np.arange(len(etot)) * 1.0

fig, axes = plt.subplots(2, 1, figsize=(12, 9), dpi=200, sharex=True)

axes[0].plot(etot, c='#3a7')
axes[0].set_ylabel('Energy [eV]')

ax1 = axes[1]
ax1.plot(time, temp, c='#c33')
ax1.set_ylabel('Temperature (kelvin)')

ax2 = ax1.twinx()
ax2.plot(time, press, c='#649')
ax2.set_ylabel('Pressure (bar)')

for ax in axes.flatten():
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='major', width=.8, length=4)
    ax.tick_params(which='minor', width=.6, length=3)
    ax.grid(color='grey', linestyle='-', linewidth=.25)

fig.suptitle(f'NVT observables for p={int(p)} bar and T={t:.1f} K')
fig.supxlabel('Time (fs)')
fig.tight_layout(pad=.8)
plt.show()
