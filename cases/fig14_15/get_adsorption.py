import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
from matplotlib import rc
import sys

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

data_dns = pd.read_csv('data_dns.csv')
data_pnm = pd.read_csv('data_pnm.csv')

coords_dns_old = data_dns['arc_length'].tolist()

left = 760
right = 840
coords_dns = coords_dns_old[:left] + coords_dns_old[right:]
coords_dns[left-1] = 0.004
coords_dns[left] = 0.0040001

coords_pnm_1 = data_pnm['X_cells'].tolist()
coords_pnm_2 = list(np.linspace(4.e-3, 5.0e-3, 10))
coords_pnm = coords_pnm_1 + coords_pnm_2
# 
pressure_dns_old = data_dns['p'].tolist()
pressure_dns = pressure_dns_old[:left] + pressure_dns_old[right:]
pressure_dns[left-1] = 3.32e+6

pressure_pnm_1 = data_pnm['P_array'].tolist()
pressure_pnm_2 = [3.32e6] * 10
pressure_pnm = pressure_pnm_1 + pressure_pnm_2
# 

srpDep_dns_old = data_dns['alpha1M'].tolist()
srpDep_dns = srpDep_dns_old[:left] + srpDep_dns_old[right:]
srpDep_dns[left-1] = 2.96


srpDep_pnm_1 = data_pnm['alpha1M'].tolist()
srpDep_pnm_2 = [0.] * 10
srpDep_pnm = srpDep_pnm_1 + srpDep_pnm_2
# 
fig_width = 3.5
y_scale = 0.9
fig, ax = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                       tight_layout=True)
# 
# ### Fake plots for legend ###
lws = [1.5, 1.5]
line_names = ['DBS', 'HybridPNM']
line_types = ['-', '--']
for i in range(len(line_names)):
    plt.plot([], [], linestyle=line_types[i], c='k',
             label=line_names[i], lw=lws[i])
legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1, bbox_to_anchor=[0., 0.8], loc='center left')
plt.gca().add_artist(legend_1)
# ### Fake plots for legend ###
# 
l1, = ax.plot(coords_dns, srpDep_dns, color='tab:orange', lw=2, alpha=0.5)
ax.plot(coords_pnm, srpDep_pnm, color='tab:orange', linestyle='--', lw=1.5)
# #
ax2 = ax.twinx()
l2, = ax2.plot(coords_dns, pressure_dns, color='tab:blue', lw=2, alpha=0.5)
ax2.plot(coords_pnm, pressure_pnm, color='tab:blue', linestyle='--', lw=1.5)
# 
ax.set_xlabel('$Y$, $m$')
ax.set_ylabel('$\\alpha$')
ax2.set_ylabel('$P$, $Pa$')
# 

line_types = ['-']
line1 = ax.plot([], [], linestyle=line_types[0], c='tab:orange', label='$\\alpha$')
line2 = ax.plot([], [], linestyle=line_types[0], c='tab:blue', label='$P$')
lns1 = line1 + line2
labs1 = [l.get_label() for l in lns1]

ax.legend(lns1, labs1, bbox_to_anchor=[0., 0.55], loc='center left')

plt.savefig('fig15_adsorption_2phase_1D.pdf', format="pdf",
            bbox_inches='tight')
plt.show()

