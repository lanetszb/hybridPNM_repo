import sys
import os
import numpy as np
import copy as cp
import math as m
from alive_progress import alive_bar, config_handler
import random
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../../'))
sys.path.append(os.path.join(current_path, '../../../netgrid'))
sys.path.append(os.path.join(current_path, '../../../raplea'))

import hybridPNM as pnm

import raplea as rl

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
# plt.style.use('dracula')

np.set_printoptions(threshold=np.inf, linewidth=300)

# ==========================================================================================
#                                       input                                              #
# ==========================================================================================

dt = 1.e-1
t_max = 1000.
t_bc = 1.e1
t_step_out_coarse = t_max / 1000
t_step_out_fine = t_max / 1000

visc0 = 1.e-3
visc1 = 2.e-5
gamma = 1
beta = 7.2e-2
theta = (m.pi / 180.) * 30
xi = 0.

P_ini = 3.0e6
PBcTypes = {"inlet": "Dirichlet", "outlet": "Dirichlet"}
valuePInlet = 4.0e6
valuePInletIni = P_ini
valuePOutlet = 4.0e6
valuePOutletIni = P_ini

S0_ini = 0.
S0BcTypes = {"inlet": "Neumann", "outlet": "Neumann"}
valueS0Inlet = 0.
valueS0Outlet = 0.

scaleM = 0.5
phiM = 0.05
kM = 4. * 1.e-19

S0RM = 0.0
S1RM = 0.0

P0M_ini = P_ini
S0M_ini = 0.

round_output_time = 4

# ==========================================================================================
#                                    netgrid for pnm                                       #
# ==========================================================================================


poresCoordinates = {0: [0., 0.005], 1: [1.e-2, 0.005]}
throatsPores = {0: [0, 1]}
throatsWidths = {0: 0.002}
throatsDepths = {0: 0.000125}
deltaV = 0.05
minCellsN = 80
inletPores = {0}
outletPores = {1}
netgrid = Netgrid(poresCoordinates, throatsPores,
                  throatsWidths, throatsDepths, deltaV, minCellsN,
                  inletPores, outletPores)

# ==========================================================================================
#                                    netgrid for matrix                                    #
# ==========================================================================================


poresCoordinatesM = {}
throatsPoresM = {}
throatsWidthsM = {}
throatsDepthsM = {}

deltaVM = scaleM * deltaV / phiM
minCellsNM = 32
inletPoresM = set()
outletPoresM = set()
throatsVM = {0: 1e-08}

poreM = 0
for throat, pores in throatsPores.items():
    width = throatsWidths[throat]
    depth = throatsDepths[throat]
    L = netgrid._throatsLs[throat]
    LM = throatsVM[throat] / depth / L * scaleM
    pore_a = pores[0]
    pore_b = pores[1]
    coordinates_a = poresCoordinates[pore_a]
    coordinates_b = poresCoordinates[pore_b]
    alpha = -m.atan2(coordinates_b[0] - coordinates_a[0], coordinates_b[1] - coordinates_a[1])
    cellsN = netgrid._throatsCellsNs[throat]
    countCell = 0
    for cell in netgrid._throatsCells[throat]:
        point_a = [coordinates_a[0] + (coordinates_b[0] - coordinates_a[0]) * countCell / cellsN,
                   coordinates_a[1] + (coordinates_b[1] - coordinates_a[1]) * countCell / cellsN]
        point_b = [coordinates_a[0] + (coordinates_b[0] - coordinates_a[0]) * (countCell + 1) / cellsN,
                   coordinates_a[1] + (coordinates_b[1] - coordinates_a[1]) * (countCell + 1) / cellsN]
        point_c = [(point_a[0] + point_b[0]) / 2. + m.cos(alpha) * width / 2.,
                   (point_a[1] + point_b[1]) / 2. + m.sin(alpha) * width / 2.]
        point_pc = [point_c[0] + m.cos(alpha) * LM,
                    point_c[1] + m.sin(alpha) * LM]

        poresCoordinatesM[poreM] = point_pc
        poresCoordinatesM[poreM + 1] = point_c
        throatsPoresM[cell] = [poreM, poreM + 1]
        throatsWidthsM[cell] = L / cellsN
        throatsDepthsM[cell] = depth
        inletPoresM.add(poreM)
        outletPoresM.add(poreM + 1)

        poreM += 2
        countCell += 1

netgridM = Netgrid(poresCoordinatesM, throatsPoresM,
                   throatsWidthsM, throatsDepthsM, deltaVM, minCellsNM,
                   inletPoresM, outletPoresM)

# ==========================================================================================
#                                        raplea                                            #
# ==========================================================================================


# alphaCoeffs = [-2.19813E-26, 3.73110E-19, 2.50417E-12, 9.35591E-06, 1.96098E-01]
alphaCoeffs = [0, 0, 0, 4.3549791928771896e-06, -10.161618116713441]

pars = rl.Pars(S0R=S0RM, S1R=S1RM, M0=1., M1=1., MC=0., m0=1., m1=1., mC=1., mu0=visc0, mu1=visc1,
               alphaCoeffs=alphaCoeffs, sorptionType='kR1Linear', is0Wet=True, dt=dt)

alpha1M_ini = pars.alpha1(P0M_ini)

throatsNM = netgridM._throatsN
facesNM = netgridM._facesN
cellsNM = netgridM._cellsN
uM = [np.zeros(facesNM, dtype=np.float64),
      np.zeros(facesNM, dtype=np.float64)]
P0M = [np.full(cellsNM, P0M_ini, dtype=np.float64),
       np.full(cellsNM, P0M_ini, dtype=np.float64)]
S0M = [np.full(cellsNM, S0M_ini, dtype=np.float64),
       np.full(cellsNM, S0M_ini, dtype=np.float64)]
alpha1M = [np.full(cellsNM, alpha1M_ini, dtype=np.float64),
           np.full(cellsNM, alpha1M_ini, dtype=np.float64)]

data = {"uM": np.zeros(throatsNM, dtype=np.float64),
        "qM": np.zeros(throatsNM, dtype=np.float64),
        "S0M": np.zeros(throatsNM, dtype=np.float64),
        "fM": np.zeros(throatsNM, dtype=np.float64),
        "deltaAlpha1M": np.zeros(throatsNM, dtype=np.float64),
        "avAlpha1M": np.full(throatsNM, alpha1M_ini, dtype=np.float64),
        "deltaS0M": np.zeros(throatsNM, dtype=np.float64),
        "avS0M": np.full(throatsNM, S0M_ini, dtype=np.float64),
        "S0": np.full(throatsNM, S0_ini, dtype=np.float64),
        "P": np.full(throatsNM, P_ini, dtype=np.float64),
        "KFactor": np.ones(throatsNM, dtype=np.float64)}

matrix = rl.Matrix(phiM, kM, scaleM, netgridM, pars, data, uM, P0M, S0M, alpha1M)

matrix.setSolverPs('sparseLU')
matrix.setSolverSs('sparseLU')

iPrevM = matrix._iPrev
uM_array_times = [np.copy(uM[iPrevM])]
uM_av_times = [np.mean(uM[iPrevM])]
P0M_array_times = [np.copy(P0M[iPrevM])]
P0M_av_times = [np.mean(P0M[iPrevM])]
S0M_array_times = [np.copy(S0M[iPrevM])]
S0M_av_times = [np.mean(S0M[iPrevM])]
alpha1M_array_times = [np.copy(alpha1M[iPrevM])]
alpha1M_av_times = [np.mean(alpha1M[iPrevM])]

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
equationP._bcTypes = PBcTypes
equationP._bc = {"inlet": valuePInlet, "outlet": valuePOutlet}
equationP._solver = 'sparseLU'

locS = pnm.LocS(netgrid, props, u, dS0, P, S0, data)
conS = pnm.ConS(netgrid, props, u, dS0, P, S0, data)
equationS = pnm.Equation(locS, conS, S0)
equationS._bcTypes = S0BcTypes
equationS._bc = {"inlet": valueS0Inlet, "outlet": valueS0Outlet}
equationS._solver = 'sparseLU'

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
P_array_times = [np.copy(P[iPrev])]
P_av_times = [np.mean(P[iPrev])]
u_array_times = [np.copy(u)]
u_av_times = [np.mean(np.fabs(u))]
dS0_array_times = [np.copy(dS0)]
dS0_av_times = [np.mean(dS0)]
S0_array_times = [np.copy(S0[iPrev])]
S0_av_times = [np.mean(S0[iPrev])]
S0S1 = np.copy(S0[iPrev]) * (1. - np.copy(S0[iPrev]))
S0S1_array_times = [np.copy(S0S1)]
S0S1_av_times = [np.mean(S0S1)]

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

        # pnm
        equationP._bc = {"inlet": smooth_bc(valuePInletIni, valuePInlet, t_bc, t_cur),
                         "outlet": smooth_bc(valuePOutletIni, valuePOutlet, t_bc, t_cur)}
        cfd_procedure()
        iPrev = equationP._iPrev
        S0S1 = np.copy(S0[iPrev]) * (1. - np.copy(S0[iPrev]))
        iPrev = equationP._iPrev

        # matrix
        matrix.imposeBc()
        matrix.solve()
        matrix.calculateData()
        matrix.iterateIndices()
        iPrevM = matrix._iPrev

        if (i + 1) % step_out_coarse == 0:
            ts_out_coarse.append(t_cur)
            # pnm
            P_array_times.append(np.copy(P[iPrev]))
            S0_array_times.append(np.copy(S0[iPrev]))
            u_array_times.append(np.copy(u))
            dS0_array_times.append(np.copy(dS0))
            S0S1_array_times.append(np.copy(S0S1))
            # matrix
            uM_array_times.append(np.copy(uM[iPrevM]))
            P0M_array_times.append(np.copy(P0M[iPrevM]))
            S0M_array_times.append(np.copy(S0M[iPrevM]))
            alpha1M_array_times.append(np.copy(alpha1M[iPrevM]))

        if (i + 1) % step_out_fine == 0:
            ts_out_fine.append(t_cur)
            # pnm
            P_av_times.append(np.mean(P[iPrev]))
            u_av_times.append(np.mean(np.fabs(u)))
            dS0_av_times.append(np.mean(dS0))
            S0_av_times.append(np.mean(S0[iPrev]))
            S0S1_av_times.append(np.mean(S0S1))
            # matrix
            uM_av_times.append(np.mean(uM[iPrevM]))
            P0M_av_times.append(np.mean(P0M[iPrevM]))
            S0M_av_times.append(np.mean(S0M[iPrevM]))
            alpha1M_av_times.append(np.mean(alpha1M[iPrevM]))

        bar()

# ==========================================================================================
#                                     pnm VTK output                                       #
# ==========================================================================================


cells_file_name = 'adsorption.pvd'
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
    cells_files_names.append(str(round(time, round_output_time)) + '.vtu')
    cells_files_descriptions.append(str(round(time, round_output_time)))
    netgrid.saveCells(cells_files_names[-1])
    saveFilesCollectionToFile(cells_file_name, cells_files_names, cells_files_descriptions)

# ==========================================================================================
#                                     matrix VTK output                                    #
# ==========================================================================================


cells_file_name = 'adsorption_M.pvd'
cells_files_names = list()
cells_files_descriptions = list()
for i in range(len(ts_out_coarse)):
    time = ts_out_coarse[i]
    netgridM.cellsArrays = {'P0M': P0M_array_times[i],
                            'S0M': S0M_array_times[i],
                            'alpha1M': alpha1M_array_times[i]}
    netgridM.facesArrays = {'uM': uM_array_times[i]}
    cells_files_names.append(str(round(time, round_output_time)) + '_M.vtu')
    cells_files_descriptions.append(str(round(time, round_output_time)))
    netgridM.saveCells(cells_files_names[-1])
    saveFilesCollectionToFile(cells_file_name, cells_files_names, cells_files_descriptions)

# ==========================================================================================
#                             matrix csv output                                            #
# ==========================================================================================    
X_cells = np.linspace(0, 4.e-3, minCellsNM)
P_array = P0M_array_times[1000][0:32]
alpha1M = alpha1M_array_times[1000][0:32]

d = {'X_cells': X_cells, 'P_array': P_array, 'alpha1M': alpha1M}
df = pd.DataFrame(data=d)
df.to_csv('data_pnm.csv')
