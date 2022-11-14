import sys
import os
import numpy as np
import copy as cp
import math as m
# from alive_progress import alive_bar, config_handler
import random
import json
import pandas as pd
from scipy import stats

import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"


# plt.style.use('dracula')

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return r_value ** 2


### openFoam output ###
data_dns = pd.read_csv('data_dns.csv')
timesDns = data_dns['timesDns']
timesMatrixAlpha1Dns = data_dns['timesMatrixAlpha1Dns']
timesSrpDepDns = data_dns['timesSrpDepDns']
timesMatrixPDns = data_dns['timesMatrixPDns']

# hybridPNM output
data_pnm = pd.read_csv('data_pnm.csv')

timesPNM = data_pnm['ts_out_fine']
timesMatrixAlpha1Pnm = data_pnm['S0M_av_times']
timesSrpDepPnm = data_pnm['alpha1M_av_times']
timesMatrixPPnm = data_pnm['P0M_av_times']

# print(len(timesMatrixAlpha1Pnm))
timesPNM.pop(0)
timesMatrixAlpha1Pnm.pop(0)
timesSrpDepPnm.pop(0)
timesMatrixPPnm.pop(0)
# R2 values
matrixAlphaR2 = rsquared(timesMatrixAlpha1Dns, timesMatrixAlpha1Pnm)
timesSrpDepR2 = rsquared(timesSrpDepDns, timesSrpDepPnm)
timesMatrixPR2 = rsquared(timesMatrixPDns, timesMatrixPPnm)

fig_width = 3.5
y_scale = 0.9
fig0, ax0 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)
# #
# # 1rd set
lws = [1.5, 1.5]
line_names = ['DBS', 'HybridPNM']
line_types = ['-', '--']
for i in range(len(line_names)):
    plt.plot([], [], linestyle=line_types[i], c='k',
             label=line_names[i], lw=lws[i])
legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1)
plt.gca().add_artist(legend_1)

ax0.plot(timesDns, timesMatrixAlpha1Dns, color='tab:red', lw=2, alpha=0.5)
ax0.plot(timesPNM, timesMatrixAlpha1Pnm, linestyle='--', color='tab:red', lw=1.5)
ax0.text(60., 0.5, '$R^{\\mathit2}=' + str(round(matrixAlphaR2, 2)) + '$')
ax0.set_ylabel('$S_{\\mathit{0}}$')
ax0.set_xlabel('$t, s$')
# ax0.legend(loc=0)
plt.tight_layout()
plt.savefig('fig19_satM.pdf', format="pdf",
            bbox_inches='tight')
plt.show()

fig1, ax1 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)
# 2nd set
for i in range(len(line_names)):
    plt.plot([], [], linestyle=line_types[i], c='k',
             label=line_names[i], lw=lws[i])
legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1)
plt.gca().add_artist(legend_1)

ax1.plot(timesDns, timesSrpDepDns, color='tab:orange', lw=2, alpha=0.5)
ax1.plot(timesPNM, timesSrpDepPnm, linestyle='--', color='tab:orange', lw=1.5)
ax1.text(60., 2.95, '$R^{\\mathit2}=' + str(round(timesSrpDepR2, 2)) + '$')
ax1.set_ylabel('$\\alpha$')
ax1.set_xlabel('$t, s$')
# ax1.legend(loc=0)
plt.tight_layout()
plt.savefig('fig20_alphaM.pdf', format="pdf",
            bbox_inches='tight')
plt.show()
#
fig2, ax2 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)

for i in range(len(line_names)):
    plt.plot([], [], linestyle=line_types[i], c='k',
             label=line_names[i], lw=lws[i])
legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc=4)
plt.gca().add_artist(legend_1)

ax2.plot(timesDns, timesMatrixPDns, color='tab:blue', lw=2, alpha=0.5)
ax2.plot(timesPNM, timesMatrixPPnm, linestyle='--', color='tab:blue', lw=1.5)
ax2.text(20., 3.8e+6, '$R^{\\mathit2}=' + str(round(timesMatrixPR2, 3)) + '$')
ax2.set_ylabel('$P, Pa$')
ax2.set_xlabel('$t, s$')

plt.tight_layout()
plt.savefig('fig21_pM.pdf', format="pdf",
            bbox_inches='tight')
plt.show()
#
