import sys
import os
import numpy as np
import copy as cp
import math as m
import pandas as pd
from alive_progress import alive_bar, config_handler
import random
import json

import matplotlib.pyplot as plt
from matplotlib import rc

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../../'))
sys.path.append(os.path.join(current_path, '../../../netgrid'))
sys.path.append(os.path.join(current_path, '../../../raplea'))

import hybridPNM as pnm

from netgrid import Netgrid, saveFilesCollectionToFile

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
plt.style.use('dracula')

np.set_printoptions(threshold=np.inf, linewidth=300)

# ==========================================================================================
#                                       input                                              #
# ==========================================================================================

u_inlet = 0.005
dt = 1.e-3
t_max = 1.
t_bc = 0.
t_step_out_coarse = t_max / 100
t_step_out_fine = t_max / 100

visc0 = 2.e-5
visc1 = 2.e-5
gamma = 1.
beta = 0.
theta = (m.pi / 180.) * 45
xi = 0.

width = 2.e-4
depth = 5.e-5

P_ini = 1.e6
PBcTypes = {"inlet": "Dirichlet", "outlet": "Dirichlet"}
valuePInlet = 0.50000001e6
valuePInletIni = 0.50000001e6
valuePOutlet = 0.49999999e6
valuePOutletIni = 0.49999999e6

S0_ini = 0.
S0BcTypes = {"inlet": "Neumann", "outlet": "Neumann"}
valueS0Inlet = 0.
valueS0Outlet = 0.

round_output_time = 7

# ==========================================================================================
#                                    netgrid for pnm                                       #
# ==========================================================================================


poresCoordinates = {0: [0., 0.], 1: [1.e-3, 0]}
throatsPores = {0: [0, 1]}
throatsWidths = {0: width}
throatsDepths = {0: depth}
deltaV = 0.5
minCellsN = 50
inletPores = {1}
outletPores = {0}
netgrid = Netgrid(poresCoordinates, throatsPores,
                  throatsWidths, throatsDepths, deltaV, minCellsN,
                  inletPores, outletPores)

# ==========================================================================================
#                                        raplea                                            #
# ==========================================================================================


throatsNM = netgrid._cellsN

data = {"uM": np.zeros(throatsNM, dtype=np.float64),
        "qM": np.zeros(throatsNM, dtype=np.float64),
        "S0M": np.zeros(throatsNM, dtype=np.float64),
        "fM": np.zeros(throatsNM, dtype=np.float64),
        "deltaAlpha1M": np.zeros(throatsNM, dtype=np.float64),
        "avAlpha1M": np.zeros(throatsNM, dtype=np.float64),
        "deltaS0M": np.zeros(throatsNM, dtype=np.float64),
        "avS0M": np.zeros(throatsNM, dtype=np.float64),
        "S0": np.zeros(throatsNM, dtype=np.float64),
        "P": np.zeros(throatsNM, dtype=np.float64),
        "KFactor": np.ones(throatsNM, dtype=np.float64)}

# ==========================================================================================
#                                     hybridPNM                                            #
# ==========================================================================================


props = pnm.Props({"dt": float(dt),
                   "visc0": float(visc0), "visc1": float(visc0),
                   "theta": float(theta), "gamma": float(gamma), "beta": float(beta),
                   "dim": "2D", "xi": float(xi)})

u = np.zeros(netgrid._facesN, dtype=np.float64)
dS0 = np.zeros(netgrid._facesN, dtype=np.float64)
P = [np.full(netgrid._cellsN, P_ini, dtype=np.float64), np.full(netgrid._cellsN, P_ini, dtype=np.float64)]
S0 = [np.full(netgrid._cellsN, S0_ini, dtype=np.float64), np.full(netgrid._cellsN, S0_ini, dtype=np.float64)]

locP = pnm.LocP(netgrid, props, u, dS0, P, S0, data)
conP = pnm.ConP(netgrid, props, u, dS0, P, S0, data)
equationP = pnm.Equation(locP, conP, P)
equationP._solver = 'sparseLU'
equationP._bcTypes = PBcTypes
equationP._bc = {"inlet": valuePInlet, "outlet": valuePOutlet}

locS = pnm.LocS(netgrid, props, u, dS0, P, S0, data)
conS = pnm.ConS(netgrid, props, u, dS0, P, S0, data)
equationS = pnm.Equation(locS, conS, S0)
equationS._solver = 'sparseLU'
equationS._bcTypes = S0BcTypes
equationS._bc = {"inlet": valueS0Inlet, "outlet": valueS0Outlet}


def cfd_procedure():
    equationP.calculateCoefficients()
    equationP.fillMatrix()
    equationP.solve()

    equationS.calculateCoefficients()
    equationS.fillMatrix()
    equationS.solve()
    equationS.calculateData()

    equationP.iterateIndices()
    equationS.iterateIndices()


iPrev = equationP._iPrev
P_array_times, P_av_times = [cp.deepcopy(P[iPrev])], [np.mean(P[iPrev])]
u_array_times, u_av_times = [cp.deepcopy(u)], [np.mean(np.fabs(u))]
dS0_array_times, dS0_av_times = [cp.deepcopy(dS0)], [np.mean(dS0)]
S0_array_times, S0_av_times = [cp.deepcopy(S0[iPrev])], [np.mean(S0[iPrev])]
S0S1 = cp.deepcopy(S0[iPrev]) * (1. - cp.deepcopy(S0[iPrev]))
S0S1_array_times, S0S1_av_times = [cp.deepcopy(S0S1)], [np.mean(S0S1)]

# ==========================================================================================
#                                     time loop                                            #
# ==========================================================================================

if t_step_out_coarse < dt:
    t_step_out_coarse = dt

if t_step_out_fine < dt:
    t_step_out_fine = dt

t_stepsN = round(t_max / dt)
step_out_coarse = round(t_step_out_coarse / dt)
step_out_fine = round(t_step_out_fine / dt)
ts_out_coarse = [0]
ts_out_fine = [0]
t_cur = 0
with alive_bar(t_stepsN, title='time loop') as bar:
    for i in range(t_stepsN):
        t_cur += dt

        pnm.modKFactor(data['KFactor'], P[iPrev])

        equationP._bc = {"inlet": smooth_bc(valuePInletIni, valuePInlet, t_bc, t_cur),
                         "outlet": smooth_bc(valuePOutletIni, valuePOutlet, t_bc, t_cur)}
        cfd_procedure()
        iPrev = equationP._iPrev
        S0S1 = cp.deepcopy(S0[iPrev]) * (1. - cp.deepcopy(S0[iPrev]))
        iPrev = equationP._iPrev

        if (i + 1) % step_out_coarse == 0:
            ts_out_coarse.append(t_cur)

            P_array_times.append(cp.deepcopy(P[iPrev]))
            S0_array_times.append(cp.deepcopy(S0[iPrev]))
            u_array_times.append(cp.deepcopy(u))
            dS0_array_times.append(cp.deepcopy(dS0))
            S0S1_array_times.append(cp.deepcopy(S0S1))

        if (i + 1) % step_out_fine == 0:
            ts_out_fine.append(t_cur)

            P_av_times.append(np.mean(P[iPrev]))
            u_av_times.append(np.mean(np.fabs(u)))
            dS0_av_times.append(np.mean(dS0))
            S0_av_times.append(np.mean(S0[iPrev]))
            S0S1_av_times.append(np.mean(S0S1))

        bar()

# ==========================================================================================
#                                     pnm VTK output                                       #
# ==========================================================================================


cells_file_name = 'compressibility.pvd'
cells_files_names = list()
cells_files_descriptions = list()
for i in range(len(ts_out_coarse)):
    time = ts_out_coarse[i]
    P = P_array_times[i]
    S0 = S0_array_times[i]
    S0S1 = S0S1_array_times[i]
    u = u_array_times[i]
    dS0 = dS0_array_times[i]
    netgrid.cellsArrays = {'P': P, 'S0': S0, 'S0S1': S0S1}
    netgrid.facesArrays = {'u': u, 'dS0': dS0}
    cells_files_names.append('compressibility_' + str(round(time, round_output_time)) + '.vtu')
    cells_files_descriptions.append(str(round(time, round_output_time)))
    netgrid.saveCells(cells_files_names[-1])
    saveFilesCollectionToFile(cells_file_name, cells_files_names, cells_files_descriptions)

u_out = u_array_times[-1][-1]
P_av = (valuePInlet + valuePOutlet) / 2.
length = poresCoordinates[1][0] - poresCoordinates[0][0]
k = (valuePInlet - valuePOutlet) * u_out * visc0 / length

data = {'P_av': P_av, 'k': k}

json_file_name = 'k_0.json'
with open(json_file_name, 'w') as f:
    json.dump(data, f, sort_keys=False, indent=4 * ' ', ensure_ascii=False)
