import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import glob
import json
import os
import sys

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

def poro_strain_stress(phi_0, P, P_0, M, S_max, P_L):
	return 1. + 1. / M / phi_0 * (P - P_0) + S_max / 3. / phi_0 * (K / M - 1) * ((P / (P_L + P)) - (P_0 / (P_L + P_0)))
    
def poro_constants(P, c_1, c_2, c_3, c_4):
    return (c_1 * P**2 + c_2 * P) / (c_3 + P) + c_4
    
# Constants physical
S_max = 0.035
P_L = 2.1e+6
E = 1.378e9
v = 0.35
K = E / 3. / (1. - 2 * v)
M = E * (1. - v) / (1. + v) / (1. - 2. * v)
# phi_0 = 0.0073
phi_0 = 0.01

# Constants parametric
c_1 = 4.5216e-8
c_2 = -0.263946
c_3 = 2.1e6
c_4 = 1.070558193548387


P_min = 0.5e6
P_max = 4.5e6
points_n = 300
P_range = np.linspace(P_min, P_max, points_n)

P_0 = 0.5e6

# poro_range = []
# for P in P_range:
#     poro_range.append(poro_constants(P, c_1, c_2, c_3, c_4))
poro_range = []
for P in P_range:
    poro_range.append(poro_strain_stress(phi_0, P, P_0, M, S_max, P_L))

	
poro_range_2 = np.array(poro_range) * np.array(poro_range)
poro_range_3 = np.array(poro_range) * np.array(poro_range) * np.array(poro_range)

fig_width = 3.5
y_scale = 0.9
fig, ax = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                       tight_layout=True)

markers = ['o', '']
lws = [0.3, 1.5]
line_names = ['numerical', 'analytical']
line_types = ['', '-']
for i in range(len(line_names)):
    plt.plot([], [], linestyle=line_types[i], c='k',
             label=line_names[i], marker=markers[i], lw=lws[i],
             markersize='2.5')
legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc=9)
plt.gca().add_artist(legend_1)

# l1, = ax.plot(P_range, poro_range, label = '$\\frac{\phi}{\phi_{\\mathit0}}$', color='tab:gray')
l2, = ax.plot(P_range, poro_range_2, label = '$\\left(\\frac{\phi}{\phi_{\\mathit0}}\\right)^{\\mathit2}$', color='tab:olive')

file_names = glob.glob('*.json')

length = 0.001
perms = dict()

for name in file_names:
    with open(name) as f:
        perm_data = json.load(f)
    P_av = perm_data['P_av']
    k = perm_data['k']

    filename_only = os.path.basename(name)
    result = {'perm': k, 'press': P_av}
    perms[filename_only] = result

ax.set_xlabel('$P$, $Pa$')
ax.set_ylabel('$\\frac{k}{k_{\\mathit0}}$')
ax.set_ylim(bottom=0.924, top=1.015)


for key, value in perms.items():
    press = value['press']
    perm = value['perm'] / perms['k_0.json']['perm']
    # perm = value['perm'] / perms['k0_3.0MPa.json']['perm']
    ax.plot(press, perm, color='tab:brown', marker='o', markersize=6, linestyle='')
    if key == list(perms.keys())[-1]:
        l3, = ax.plot(press, perm, color='tab:brown', marker='o', markersize=6, linestyle='', label='$\\frac{k}{k_{\\mathit0}}$')
    
# ax.legend(handles=[l1, l2, l3], bbox_to_anchor=(0.45, 0.5725), loc=3)

plt.savefig('fig9_compressibility.pdf', format="pdf",
            bbox_inches='tight')
plt.show()