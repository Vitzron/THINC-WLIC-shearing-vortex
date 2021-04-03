import numpy as np
import matplotlib.pyplot as plt
import sys

X       = np.loadtxt('X.dat')
Y       = np.loadtxt('Y.dat')

phi_files = ['phi_0.dat', 'phi_2.dat', 'phi_4.dat', 'phi_6.dat', 'phi_8.dat']
t_data = [0.0, 2.0, 4.0, 6.0, 8.0]
print(phi_files[0])

levels = np.linspace(0, 1, 11)
fig, axes = plt.subplots(1, 5, figsize=(15, 3))
for i in np.arange(len(phi_files)):
    alpha = np.loadtxt(phi_files[i])
    alpha = np.where(alpha > 0, alpha, 0)
    alpha = np.where(alpha < 1, alpha, 1)
    axes[i].set_xlim(0, 1)
    axes[i].set_ylim(0, 1)
    axes[i].set_xticks([])
    axes[i].set_yticks([])
    t = '%.1f' % t_data[i]
    axes[i].set_title('$t=$' + t + r'$\,\mathrm{s}$', fontsize='12')
    axes[i].annotate('$t=$' + t + r'$\,\mathrm{s}$', xy=(0.5,0.05), fontsize=12, va='center', ha='center')
    axes[i].contourf(X, Y, alpha, levels, cmap=plt.cm.bwr)
    plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0, hspace=0, wspace=0)

if len(sys.argv) > 1:
    plt.savefig(sys.argv[1], format='png', dpi=300)

plt.show()