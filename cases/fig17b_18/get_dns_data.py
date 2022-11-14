import sys
import os
import numpy as np
import copy as cp
import math as m
from alive_progress import alive_bar, config_handler
import random
import json

import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"


with open('data_dns.json') as f:
    data_dns = json.load(f)
    
times = sorted(data_dns['timesAlpha1'].keys(), key=lambda x: float(x[0]))
times = [float(x) for x in times]

timesMatrixAlpha1 = data_dns['timesMatrixAlpha1'].items()
timesMatrixAlpha1 = [seq[1] for seq in timesMatrixAlpha1]

timesVoidAlpha1 = data_dns['timesVoidAlpha1'].items()
timesVoidAlpha1 = [seq[1] for seq in timesVoidAlpha1]

timesSrpDep = data_dns['timesSrpDep'].items()
timesSrpDep = [seq[1] for seq in timesSrpDep]
timesSrpDep = np.array(timesSrpDep) / 20.666
timesSrpDep = list(timesSrpDep)

timesP = data_dns['timesP'].items()
timesP = [seq[1] for seq in timesP]

timesMatrixP = data_dns['timesMatrixP'].items()
timesMatrixP = [seq[1] for seq in timesMatrixP]

timesVoidP = data_dns['timesVoidP'].items()
timesVoidP = [seq[1] for seq in timesVoidP]

print(type(timesMatrixAlpha1))
        
fig_width = 12
fig, axs = plt.subplots(3, 2, figsize=(fig_width, fig_width))
# 
# # 1th set
ax0 = axs[0, 0]
ax0.plot(times, timesVoidAlpha1, label='$S_{\\mathit{0}av}$', marker='o', markersize=2, color='tab:green')
ax0.set_ylabel('$S_{\\mathit{0}}$')
ax0.set_xlabel('$t, s$')
ax0.legend(loc=0)
# 
# # 2th set
ax2 = axs[1, 0]
ax2.plot(times, timesP, label='$P_{av}$', marker='o', markersize=2, color='tab:cyan')
ax2.set_ylabel('$P, Pa$')
ax2.set_xlabel('$t, s$')
ax2.legend(loc=0)
# 
# 3rd set
ax4 = axs[2, 0]
ax4.plot(times, timesMatrixAlpha1, label='$S_{\\mathit{0}avM}$', marker='o', markersize=2, color='tab:green')
ax4.set_ylabel('$S_{\\mathit{0}}$')
ax4.set_xlabel('$t, s$')
ax4.legend(loc=0)

# 6th set
ax5 = axs[0, 1]
ax5.plot(times, timesSrpDep, label='$\\alpha_{\\mathit{1}avM}$', marker='o', markersize=2, color='tab:red')
ax5.set_ylabel('$\\alpha_{\\mathit{1}}$')
ax5.set_xlabel('$t, s$')
ax5.legend(loc=0)

# 7th set
ax6 = axs[1, 1]
ax6.plot(times, timesMatrixP, label='$P_{\\mathit{0}avM}$', marker='o', markersize=2, color='tab:cyan')
ax6.set_ylabel('$P_{\\mathit{0}M}, Pa$')
ax6.set_xlabel('$t, s$')
ax6.legend(loc=0)

# 8th set
ax7 = axs[2, 1]
ax7.plot(times, timesVoidP, label='$P_{\\mathit{0}avV}$', marker='o', markersize=2, color='tab:orange')
ax7.set_ylabel('$P_{\\mathit{0}V}, Pa$')
ax7.set_xlabel('$t, s$')
ax7.legend(loc=0)

plt.tight_layout()
plt.show()
