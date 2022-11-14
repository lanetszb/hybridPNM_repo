from pandas import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import json
import sys
import os

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

data_vof = read_csv('capillary_vof_imb.csv')
x_vof = data_vof['arc_length'].tolist()
sat_vof = data_vof['alpha.water'].tolist()
p_vof = data_vof['p_rgh'].tolist()

data_pnm = read_csv('capillary_pnm_imb.csv')
x_pnm = data_pnm['X_cells'].tolist()
sat_pnm = data_pnm['S0_array'].tolist()
p_pnm = data_pnm['P_array'].tolist()

fig_width = 3.5
y_scale = 0.9
fig, ax = plt.subplots(figsize=(fig_width, fig_width * y_scale))

line_names = ['VOF', 'HybridPNM']
line_types = ['-', '--']
line1, = ax.plot([], [], linestyle=line_types[0], c='k', label=line_names[0])
line2, = ax.plot([], [], linestyle=line_types[1], c='k', label=line_names[1])
legend_1 = ax.legend(handles=[line1, line2], scatterpoints=1, frameon=True, labelspacing=1, bbox_to_anchor=(1.012, 0.965))
ax.add_artist(legend_1)

lns02 = ax.plot(x_pnm, p_pnm, linestyle='--', lw=1.5, color='tab:blue', dashes=(5, 5))
lns01 = ax.plot(x_vof, p_vof, label='$P$', color='tab:blue', lw=2, alpha=0.5)

ax.set_xlabel("$X, m$")
ax.set_ylabel("$P, Pa$")
ax.set_ylim(bottom=80, top=300)
ax.set_xlim(0, 0.005)

twin1 = ax.twinx()

twin1.set_ylabel("$S_{\\mathit0}$")

lns12 = twin1.plot(x_pnm, sat_pnm, linestyle='--', lw=1.5, color='tab:red', dashes=(5, 5))
lns11 = twin1.plot(x_vof, sat_vof, label='$S_{\\mathit0}$', color='tab:red', lw=2, alpha=0.5)

line_types = ['-']
line1 = ax.plot([], [], linestyle=line_types[0], c='tab:blue', label='$P$')
line2 = ax.plot([], [], linestyle=line_types[0], c='tab:red', label='$S_{\\mathit0}$')
lns1 = line1 + line2
labs1 = [l.get_label() for l in lns1]
ax.legend(lns1, labs1, bbox_to_anchor=(0.18, 0.965))

plt.savefig('fig5_capillary_imb_1d.pdf', format="pdf",
            bbox_inches='tight')

plt.show()