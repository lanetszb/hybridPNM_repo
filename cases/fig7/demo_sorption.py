import numpy as np
import sys
import os
import copy as cp
from alive_progress import alive_bar, config_handler
import random

import matplotlib.pyplot as plt
from matplotlib import rc

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../../'))
sys.path.append(os.path.join(current_path, '../../../netgrid'))
sys.path.append(os.path.join(current_path, '../../../raplea'))

import raplea as rl
from plot_pars import plot_pars


bar_type = random.choice(['smooth', 'classic', 'brackets', 'blocks', 'bubbles', 'solid', 'checks',
                          'circles', 'squares', 'ruler', 'ruler2', 'fish', 'scuba'])
spinner_type = random.choice(['classic', 'stars', 'twirl', 'twirls', 'horizontal', 'vertical',
                              'waves', 'waves2', 'waves3', 'dots', 'triangles', 'brackets',
                              'arrow', 'arrows', 'arrows2', 'arrows_in', 'arrows_out', 'radioactive'])
config_handler.set_global(bar=bar_type, spinner=spinner_type)


def smooth_bc(x_ini, x_bc, time, time_curr):
    if time_curr > time:
        return x_bc
    else:
        return x_ini + (x_bc - x_ini) * time_curr / time


rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
# plt.style.use('dracula')

L = 0.1
A = 0.05
N = 100

S0R = 0.1
S1R = 0.1

# desorption
S0_ini = 0.9
P0_ini = 1.e6
# S0_ini = 0.1
# P0_ini = 1.e6

S0_right = 1.
# P0_right = P0_ini * 1.03
# desorption
P0_right = P0_ini * 0.97

dt = 0.001
t_max = 2.8
t_bc = 0.03
t_step_out_coarse = t_max / 7
t_step_out_fine = t_max / 1000

# alphaCoeffs = [0, 0, 0, 0, 0]
alphaCoeffs = [-2.19813E-26, 3.73110E-19, 2.50417E-12, 9.35591E-06, 1.96098E-01]

grid = rl.Grid(phi=0.2, k=1E-12, L=L, A=A, N=N)

pars = rl.Pars(S0R=S0R, S1R=S1R, M0=0.465, M1=0.86, MC=0., m0=1.709, m1=2., mC=1., mu0=2.E-5, mu1=5.E-4,
               alphaCoeffs=alphaCoeffs, sorptionType='kR1Linear', is0Wet=True, dt=dt)

# S0 = np.linspace(0, 1, 300, dtype=float)
# plot_pars(pars, S0)
# sys.exit(0)

alpha1_ini = pars.alpha1(P0_ini)

u = [np.zeros(N + 1, dtype=np.float64), np.zeros(N + 1, dtype=np.float64)]
P0 = [np.full(N, P0_ini, dtype=np.float64), np.full(N, P0_ini, dtype=np.float64)]
S0 = [np.full(N, S0_ini, dtype=np.float64), np.full(N, S0_ini, dtype=np.float64)]
alpha1 = [np.full(N, alpha1_ini, dtype=np.float64), np.full(N, alpha1_ini, dtype=np.float64)]

locP = rl.LocP(grid, pars, u, P0, S0, alpha1)
conP = rl.ConP(grid, pars, u, P0, S0, alpha1)
equationP = rl.Equation(locP, conP, X=P0)
equationP._solver = "sparseLU"
# equationP._tolerance = 1.e-20


locS = rl.LocS(grid, pars, u, P0, S0, alpha1)
conS = rl.ConS(grid, pars, u, P0, S0, alpha1)
equationS = rl.Equation(locS, conS, X=S0)
equationS._solver = "sparseLU"
# equationS._tolerance = 1.e-20
equationS._bcTypes = {'left': 'Neumann', 'right': 'Neumann'}
equationS._bc = {'left': 0., 'right': S0_right}


def cfd_sorption(equationP, equationS, P0_right, S0_right):
    equationP._bcTypes = {'left': 'Neumann', 'right': 'Dirichlet'}
    equationP._bc = {'left': 0, 'right': P0_right}
    equationS._bcTypes = {'left': 'Neumann', 'right': 'Neumann'}
    equationS._bc = {'left': 0., 'right': S0_right}

    equationP.calculateCoefficients()
    equationP.fillMatrix()
    equationP.solve()
    equationP.calculateAlpha1()
    equationP.calculateU()
    equationS.calculateCoefficients()
    equationS.fillMatrix()
    equationS.solve()
    equationP.iterateIndices()
    equationS.iterateIndices()


iPrev = equationP._iPrev
S0_array_times = [np.copy(S0[iPrev])]
S0_av_times = [np.mean(S0[iPrev])]
alpha1_array_times = [np.copy(alpha1[iPrev])]
alpha1_av_times = [np.mean(alpha1[iPrev])]
P0_array_times = [np.copy(P0[iPrev])]
P0_av_times = [np.mean(P0[iPrev])]
u_array_times = [np.copy(u[iPrev])]
u_av_times = [np.mean(u[iPrev])]

if t_step_out_coarse < dt:
    t_step_out_coarse = dt

if t_step_out_fine < dt:
    t_step_out_fine = dt

t_stepsN = round(t_max / dt)
step_out_coarse = round(t_step_out_coarse / dt)
step_out_fine = round(t_step_out_fine / dt)
ts_out_coarse = [0]
ts_out_fine = [0]
t_cur = 0.
with alive_bar(t_stepsN, title='time loop') as bar:
    for i in range(t_stepsN):
        t_cur += dt
        cfd_sorption(equationP, equationS, smooth_bc(P0_ini, P0_right, t_bc, t_cur), S0_right)

        iPrev = equationP._iPrev

        if (i + 1) % step_out_coarse == 0:
            ts_out_coarse.append(t_cur)
            S0_array_times.append(np.copy(S0[iPrev]))
            alpha1_array_times.append(np.copy(alpha1[iPrev]))
            P0_array_times.append(np.copy(P0[iPrev]))
            u_array_times.append(np.copy(u[iPrev]))

        if (i + 1) % step_out_fine == 0:
            ts_out_fine.append(t_cur)
            S0_av_times.append(np.mean(S0[iPrev]))
            alpha1_av_times.append(np.mean(alpha1[iPrev]))
            P0_av_times.append(np.mean(P0[iPrev]))
            u_av_times.append(np.mean(u[iPrev]))

        bar()

fig_width = 3.5
y_scale = 2.7
fig, axs = plt.subplots(3, 1, figsize=(fig_width, fig_width * y_scale))

X_cells = np.linspace(0, L, N)
X_faces = np.linspace(-grid.dL() / 2., L + grid.dL() / 2., N + 1)

time_round = 7

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
          'tab:purple', 'tab:gray', 'tab:pink', 'tab:olive',
          'tab:cyan']

# 1st set
ax0 = axs[0]
for i in range(len(alpha1_array_times)):
    ax0.plot(X_cells, S0_array_times[i], color=colors[i], label=str(round(ts_out_coarse[i], time_round)))
ax0.set_xlabel('$X, m$')
ax0.set_ylabel('$S_{\\mathit0}$')
ax0.set_xlim(0., L)
# ax0.set_ylim(0., 1.)
# ax0.legend(title='time, $s$')

fig.legend(title='$t, s$', bbox_to_anchor=(0.11, 0.84, 0.8, 0.15), loc=2, ncol=3, mode="expand", borderaxespad=0.)

# 2nd set
ax1 = axs[1]
for i in range(len(alpha1_array_times)):
    ax1.plot(X_cells, alpha1_array_times[i], color=colors[i], label=str(round(ts_out_coarse[i], time_round)))
ax1.set_xlabel('$X, m$')
ax1.set_ylabel('$\\alpha_\\mathit{1}$')
ax1.set_xlim(0., L)
# ax1.legend(title='time, $s$')

# 3nd set
ax2 = axs[2]
for i in range(len(alpha1_array_times)):
    ax2.plot(X_cells, P0_array_times[i], color=colors[i], label=str(round(ts_out_coarse[i], time_round)))
ax2.set_ylabel('$P, Pa$')
ax2.set_xlabel('$X, m$')
ax2.set_xlim(0., L)

plt.tight_layout()
plt.subplots_adjust(left=0.1)
fig.subplots_adjust(top=0.895)

plt.savefig('fig7_desorption.pdf', format="pdf",
            bbox_inches='tight')
plt.show()
