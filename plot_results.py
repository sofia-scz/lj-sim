import numpy as np
from process import get_everything
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   aux routines
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def sample_averages(x, ncycles, sample_size):

    ndata = len(x)
    batch_size = ndata // sample_size
    xm, x2m = np.zeros((ncycles, batch_size)), np.zeros((ncycles, batch_size))

    for i in range(ncycles):
        np.random.shuffle(x)
        for j in range(batch_size):
            xm[i, j] = x[j*sample_size:(j+1)*sample_size].mean()
            x2m[i, j] = (x[j*sample_size:(j+1)*sample_size]**2).mean()

    return xm, x2m

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    read data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


data_dir = 'results'
folders = os.listdir(data_dir)

npoints = len(folders)
data = np.zeros((npoints, 17))
npart = 10976

# units
kB = 8.617333262e-2
bar_meVangs = 1.6021766341182764e3
meVpart_Jmol = 9.3648533213043279e1
cycles = 30
sample_size = 50

for i, folder in enumerate(folders):
    temp = float(folder[:-1])
    energy, boxlen, density, pressure, _, _ = np.loadtxt(
        os.path.join(data_dir, folder, 'table.out'), skiprows=1).transpose()

    volume = boxlen**3
    pm = pressure.mean()/bar_meVangs
    enthalpy = energy + (volume/npart) * pm

    vm, v2m = sample_averages(volume, cycles, sample_size)
    bulk = vm / (v2m - vm**2) * (kB*temp) * (bar_meVangs * 1e-3)
    bulk_means, bulk_stds = bulk.mean(axis=1), bulk.std(axis=1)

    um, u2m = sample_averages(energy, cycles, sample_size)
    cv = (u2m - um**2) / (kB*temp**2) * meVpart_Jmol * npart
    cv_means, cv_stds = cv.mean(axis=1), cv.std(axis=1)

    hm, h2m = sample_averages(enthalpy, cycles, sample_size)
    cp = (h2m - hm**2) / (kB*temp**2) * meVpart_Jmol * npart
    cp_means, cp_stds = cp.mean(axis=1), cp.std(axis=1)

    alpha = (pm * (v2m - vm**2) + (um - um.mean())*(vm - vm.mean())) / vm / (kB*temp**2)
    alpha_means, alpha_stds = alpha.mean(axis=1), alpha.std(axis=1)

    data[i, 0] = temp
    data[i, 1] = energy.mean()
    data[i, 2] = energy.std()
    data[i, 3] = enthalpy.mean()
    data[i, 4] = enthalpy.std()
    data[i, 5] = density.mean()
    data[i, 6] = density.std()
    data[i, 7] = bulk_means.mean()
    data[i, 8] = bulk_stds.mean()
    data[i, 9] = alpha_means.mean()
    data[i, 10] = alpha_stds.mean()
    data[i, 11] = cp_means.mean()
    data[i, 12] = cp_stds.mean()
    data[i, 13] = cv_means.mean()
    data[i, 14] = cv_stds.mean()
    data[i, 15] = volume.mean()
    data[i, 16] = volume.std()

data = data[data[:, 0].argsort()]
x = data[:, 0]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 plot varianzas
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# fig, ax = plt.subplots(1, 1, figsize=(7, 5), dpi=200, sharex=True)

# # plot energy
# y = data[:, 2]**2
# color = '#396'
# ax.scatter(x, y, marker='D', color='#fff', edgecolor=color, zorder=1, label='u')
# ax.plot(x, y, color=color, lw=0.7, zorder=2)

# y = data[:, 4]**2
# color = '#369'
# ax.scatter(x, y, marker='D', color='#fff', edgecolor=color, zorder=1, label='h a p fijo')
# ax.plot(x, y, color=color, lw=0.7, zorder=2)

# y = data[:, -1]**2
# color = '#c33'
# ax.scatter(x, y, marker='D', color='#fff', edgecolor=color, zorder=1, label='h a p instantáneo')
# ax.plot(x, y, color=color, lw=0.7, zorder=2)
# ax.set_ylabel('Varianza (meV)')
# ax.legend()


# ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.tick_params(which='both', direction='in', top=True, right=True)
# ax.tick_params(which='major', width=.8, length=4)
# ax.tick_params(which='minor', width=.6, length=3)
# ax.grid(color='grey', linestyle='-', linewidth=.25)

# fig.supxlabel('Temperatura (K)')
# fig.tight_layout(pad=.8, w_pad=0.3, h_pad=0.2)
# plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                      plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig, axes = plt.subplots(3, 2, figsize=(12, 10), dpi=200, sharex=True)

# plot energy
ax = axes[0, 0]
y = data[:, 1]
color = '#396'
ax.scatter(x, y, marker='D', color='#fff', edgecolor=color, zorder=1)
ax.plot(x, y, color=color, lw=0.7, zorder=2)
ax.set_ylabel('Energía por átomo (meV)')

# plot cv
ax = axes[0, 1]
y = data[:, 9]
yerr = data[:, 10]
color = '#963'
skip = 1
ax.scatter(x[skip:], y[skip:], marker='s', color='#fff', edgecolor=color, zorder=1)
ax.errorbar(x[skip:], y[skip:], yerr=yerr[skip:], color=color, lw=0.7, capsize=1.5, zorder=2)
ax.set_ylabel('Coeficiente de expansión térmica (1 / K)')

# plot enthalpy
ax = axes[1, 0]
y = data[:, 3]
color = '#369'
ax.scatter(x, y, marker='^', color='#fff', edgecolor=color, zorder=1)
ax.plot(x, y, color=color, lw=0.7, zorder=2)
ax.set_ylabel('Entalpía por átomo (meV)')

# plot cp
ax = axes[1, 1]
y = data[:, 11]
yerr = data[:, 12]
color = '#993'
ax.scatter(x, y, marker='s', color='#fff', edgecolor=color, zorder=1)
ax.errorbar(x, y, yerr=yerr, color=color, lw=0.7, capsize=1.5, zorder=2)
ax.set_ylabel('Calor específico isobárico (J / K mol)')

# plot density
ax = axes[2, 0]
y = data[:, 5]
color = '#c33'
ax.scatter(x, y, marker='o', color='#fff', edgecolor=color, zorder=1)
ax.plot(x, y, color=color, lw=0.7, zorder=2)
ax.set_ylabel('Densidad (Kg/lt)')

# plot bulk
ax = axes[2, 1]
y = data[:, 7]
yerr = data[:, 8]
color = '#b79'
ax.scatter(x, y, marker='s', color='#fff', edgecolor=color, zorder=1)
ax.errorbar(x, y, yerr=yerr, color=color, lw=0.7, capsize=1.5, zorder=2)
ax.set_ylabel('Módulo de bulk (kbar)')

for i, ax in enumerate(axes.flatten()):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='major', width=.8, length=4)
    ax.tick_params(which='minor', width=.6, length=3)
    ax.grid(color='grey', linestyle='-', linewidth=.25)
    if i % 2 == 1:
        ax.tick_params(axis='y', which='both', labelleft=False, labelright=True)
        ax.yaxis.set_label_position("right")

# fig.suptitle(f'NpT observables for p={int(p)} bar and T={t:.1f} K')
fig.supxlabel('Temperatura (K)')
fig.tight_layout(pad=.8, w_pad=0.3, h_pad=0.0)
plt.show()
