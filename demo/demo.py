import sys
import os
import numpy as np
import copy as cp
import math as m
from alive_progress import alive_bar, config_handler
import random

import matplotlib.pyplot as plt
from matplotlib import rc

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))
sys.path.append(os.path.join(current_path, '../../netgrid'))
sys.path.append(os.path.join(current_path, '../../raplea'))

import hybridPNM as pnm

# raplea can be found https://github.com/AleksZhuravlyov/raplea
import raplea as rl

# netgrid can be found https://github.com/AleksZhuravlyov/raplea
from netgrid import Netgrid, saveFilesCollectionToFile

bar_type = random.choice(['smooth', 'classic', 'brackets', 'blocks', 'bubbles', 'solid', 'checks',
                          'circles', 'squares', 'ruler', 'ruler2', 'fish', 'scuba'])
spinner_type = random.choice(['classic', 'stars', 'twirl', 'twirls', 'horizontal', 'vertical',
                              'waves', 'waves2', 'waves3', 'dots', 'triangles', 'brackets',
                              'arrow', 'arrows', 'arrows2', 'arrows_in', 'arrows_out', 'radioactive'])
config_handler.set_global(bar=bar_type, spinner=spinner_type)


# smooth boundary conditions linearly regarding time
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

dt = 0.01  # computational time step
t_max = 15.  # maximum time
t_bc = 2.  # time value for boundary smoth
t_step_out_coarse = t_max / 150  # coarse time step for data output
t_step_out_fine = t_max / 300  # fine time step for data output

visc0 = 1  # viscosity of phase 0
visc1 = 1  # viscosity of phase 1
gamma = 1  # tuning parameter (exponent) for capillary pressure
beta = 1.e2  # interfacial tension
theta = (m.pi / 180.) * 60  # contact angle
xi = 0.  # tuning parameter for interface compression

P_ini = 1.e6  # initial pressure
PBcTypes = {"inlet": "Neumann", "outlet": "Dirichlet"}  # pressure boundary conditions types
valuePInlet = 1.  # value of inlet pressure boundary condition
valuePInletIni = 0.  # initial value (for smooth boundary conditions) of inlet pressure boundary condition
valuePOutlet = P_ini  # value of outlet pressure boundary condition
valuePOutletIni = P_ini  # initial value (for smooth boundary conditions) of outlet pressure boundary condition

S0_ini = 0.  # initial saturation of phase 0
S0BcTypes = {"inlet": "Neumann", "outlet": "Neumann"}  # phase 0 saturation boundary conditions types
valueS0Inlet = 1.  # value of inlet phase 0 saturation boundary condition
valueS0Outlet = 1.  # value of outlet sphase 0 aturation boundary condition

scaleM = 4.  # visualisation scale for matrix
phiM = 0.2  # matrix porosity
kM = 1.e-12  # matrix permeability

S0RM = 0.1  # irreducible phase 0 saturation for matrix
S1RM = 0.1  # irreducible phase 1 saturation for matrix

P0M_ini = P_ini  # initial matrix pressure
S0M_ini = 0.  # initial phase 0 saturation

round_output_time = 4  # number of decimal places to round time

# ==========================================================================================
#                                    netgrid for pnm                                       #
# ==========================================================================================


poresCoordinates = {0: [0., 0.], 1: [2., 0.], 2: [3.5, 1.], 3: [4., -1.]}  # pores (nodes) coordinates
throatsPores = {0: [0, 1], 1: [1, 2], 2: [1, 3]}  # couples of pores to construct throats (channels)
throatsWidths = {0: 0.1, 1: 0.15, 2: 0.25}  # throats (channels) width
throatsDepths = {0: 0.45, 1: 0.35, 2: 0.6}  # throats (channels) depth
deltaV = 0.005  # cell volume
minCellsN = 10  # minimum number of cells per throat (channel)
inletPores = {0}  # inlet type pores
outletPores = {2, 3}  # outlet type pores

# mesh
netgrid = Netgrid(poresCoordinates, throatsPores,
                  throatsWidths, throatsDepths, deltaV, minCellsN,
                  inletPores, outletPores)

# ==========================================================================================
#                                    netgrid for matrix                                    #
# ==========================================================================================


poresCoordinatesM = {}  # pores (nodes) coordinates for matrix
throatsPoresM = {}  # couples of pores to construct throats (channels) for matrix
throatsWidthsM = {}  # throats (channels) width for matrix
throatsDepthsM = {}  # throats (channels) depth for matrix

deltaVM = scaleM * deltaV / phiM  # cell volume for matrix
minCellsNM = 100  # minimum number of cells per throat (channel) for matrix
inletPoresM = set()  # inlet type pores for matrix
outletPoresM = set()  # outlet type pores for matrix
throatsVM = {0: 0.11, 1: 0.07, 2: 0.25}  # volumes per throat (channel) for matrix

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

# mesh for matrix
netgridM = Netgrid(poresCoordinatesM, throatsPoresM,
                   throatsWidthsM, throatsDepthsM, deltaVM, minCellsNM,
                   inletPoresM, outletPoresM)

# ==========================================================================================
#                                        raplea                                            #
# ==========================================================================================

# polynomial coefficients for langmuir isotherm
alphaCoeffs = [-2.19813E-26, 3.73110E-19, 2.50417E-12, 9.35591E-06, 1.96098E-01]
# alphaCoeffs = [0, 0, 0, 0, 0]

# different parameters
# S0R is irreducible phase 0 saturation
# S1R is irreducible phase 1 saturation
# M0 is magnitude for phase 0 relative permeability corry correlation
# M1 is magnitude for phase 1 relative permeability corry correlation
# MC is magnitude for capillary pressure relative permeability corry correlation
# m0 is exponent for phase 0 relative permeability corry correlation
# m1 is exponent for phase 1 relative permeability corry correlation
# mC is exponent for capillary pressure corry correlation
# mu0 is phase 0 viscosity
# mu1 is phase 1 viscosity
# alphaCoeffs is polynomial coefficients for langmuir isotherm
# sorptionType is type of sorption account ("f", "kR1", and "kR1Linear")
# is0Wet is bool flag phase 0 is wet
# dt is computational time step
pars = rl.Pars(S0R=S0RM, S1R=S1RM, M0=0.465, M1=0.86, MC=10., m0=1.709, m1=2., mC=2., mu0=2.E-5, mu1=5.E-4,
               alphaCoeffs=alphaCoeffs, sorptionType='kR1Linear', is0Wet=True, dt=dt)

alpha1M_ini = pars.alpha1(P0M_ini)  # initial phase 1 adsorbed concentration for matrix

throatsNM = netgridM._throatsN  # number of throats (channels) for matrix
facesNM = netgridM._facesN  # number of faces for matrix
cellsNM = netgridM._cellsN  # number of cells for matrix
# matrix velocity
uM = [np.zeros(facesNM, dtype=np.float64),
      np.zeros(facesNM, dtype=np.float64)]
# matrix phase 0 pressure
P0M = [np.full(cellsNM, P0M_ini, dtype=np.float64),
       np.full(cellsNM, P0M_ini, dtype=np.float64)]
# matrix phase 0 saturation
S0M = [np.full(cellsNM, S0M_ini, dtype=np.float64),
       np.full(cellsNM, S0M_ini, dtype=np.float64)]
# matrix phase 1 adsorbed concentration
alpha1M = [np.full(cellsNM, alpha1M_ini, dtype=np.float64),
           np.full(cellsNM, alpha1M_ini, dtype=np.float64)]

# structure for exchange between submodules
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

# mesh for matrix
matrix = rl.Matrix(phiM, kM, scaleM, netgridM, pars, data, uM, P0M, S0M, alpha1M)

iPrevM = matrix._iPrev  # previous time index for matrix
uM_array_times = [np.copy(uM[iPrevM])]  # velocity time array for matrix
uM_av_times = [np.mean(uM[iPrevM])]  # average velocity time array for matrix
P0M_array_times = [np.copy(P0M[iPrevM])]  # phase 0 pressure time array for matrix
P0M_av_times = [np.mean(P0M[iPrevM])]  # phase 0 average pressure time array for matrix
S0M_array_times = [np.copy(S0M[iPrevM])]  # phase 0 saturation time array for matrix
S0M_av_times = [np.mean(S0M[iPrevM])]  # phase 0 average saturation time array for matrix
alpha1M_array_times = [np.copy(alpha1M[iPrevM])]  # adsorbed phase 1 concentration time array for matrix
alpha1M_av_times = [np.mean(alpha1M[iPrevM])]  # adsorbed phase 1 average concentration time array for matrix

# ==========================================================================================
#                                     hybridPNM                                            #
# ==========================================================================================

# different properties
# dt is computational time step
# visc0 is viscosity of phase 0
# visc1 is viscosity of phase 1
# theta is contact angle
# gamma is tuning parameter (exponent) for capillary pressure
# beta is interfacial tension
# dim is dimension ("2D" or "3D")
# xi is tuning parameter for interface compression
props = pnm.Props({"dt": float(dt),
                   "visc0": float(visc0), "visc1": float(visc0),
                   "theta": float(theta), "gamma": float(gamma), "beta": float(beta),
                   "dim": "2D", "xi": float(xi)})

# velocity
u = np.zeros(netgrid._facesN, dtype=np.float64)
# delta phase 0 saturation with respect to a face
dS0 = np.zeros(netgrid._facesN, dtype=np.float64)
# pressure
P = [np.full(netgrid._cellsN, P_ini, dtype=np.float64), np.full(netgrid._cellsN, P_ini, dtype=np.float64)]
# saturation of phase 0
S0 = [np.full(netgrid._cellsN, S0_ini, dtype=np.float64), np.full(netgrid._cellsN, S0_ini, dtype=np.float64)]

locP = pnm.LocP(netgrid, props, u, dS0, P, S0, data)  # location (volume) terms for pressure equation
conP = pnm.ConP(netgrid, props, u, dS0, P, S0, data)  # convection (surface) terms for pressure equation
equationP = pnm.Equation(locP, conP, P)  # pressure equation
equationP._bcTypes = PBcTypes  # boundary conditions types for pressure equation
equationP._bc = {"inlet": valuePInlet, "outlet": valuePOutlet}  # boundary condition values for pressure equation

locS = pnm.LocS(netgrid, props, u, dS0, P, S0, data)  # location (volume) terms for phase 0 saturation equation
conS = pnm.ConS(netgrid, props, u, dS0, P, S0, data)  # convection (surface) terms for phase 0 saturation equation
equationS = pnm.Equation(locS, conS, S0)  # phase 0 saturation equation
equationS._bcTypes = S0BcTypes  # boundary conditions types for phase 0 saturation equation
equationS._bc = {"inlet": valueS0Inlet,
                 "outlet": valueS0Outlet}  # boundary condition values for phase 0 saturation equation

# procedure for one time step numerical simulation
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


cells_file_name = 'demo.pvd'
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


cells_file_name = 'demo_M.pvd'
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
#                             pnm and matrix matplotlib output                             #
# ==========================================================================================


fig_width = 12
fig, axs = plt.subplots(4, 2, figsize=(fig_width, fig_width))

# 1th set
ax0 = axs[0, 0]
ax0.plot(ts_out_fine, S0_av_times, label='$S_{\\mathit{0}av}$', marker='o', markersize=2, color='tab:green')
ax0.set_ylabel('$S_{\\mathit{0}}$')
ax0.set_xlabel('$t, s$')
ax0.legend(loc=0)

# 2th set
ax1 = axs[1, 0]
ax1.plot(ts_out_fine, dS0_av_times, label='$dS_{\\mathit{0}av}$', marker='o', markersize=2, color='tab:blue')
ax1.set_ylabel('$dS_{\\mathit{0}av}$')
ax1.set_xlabel('$t, s$')
ax1.legend(loc=0)

# 3th set
ax2 = axs[2, 0]
ax2.plot(ts_out_fine, P_av_times, label='$P_{av}$', marker='o', markersize=2, color='tab:cyan')
ax2.set_ylabel('$P, Pa$')
ax2.set_xlabel('$t, s$')
ax2.legend(loc=0)

# 4th set
ax3 = axs[3, 0]
ax3.plot(ts_out_fine, u_av_times, label='$u_{av}$', marker='o', markersize=2, color='tab:orange')
ax3.set_ylabel('$u, m/s$')
ax3.set_xlabel('$t, s$')
ax3.legend(loc=0)

# 5th set
ax4 = axs[0, 1]
ax4.plot(ts_out_fine, S0M_av_times, label='$S_{\\mathit{0}avM}$', marker='o', markersize=2, color='tab:green')
ax4.set_ylabel('$S_{\\mathit{0}}$')
ax4.set_xlabel('$t, s$')
ax4.legend(loc=0)

# 6th set
ax5 = axs[1, 1]
ax5.plot(ts_out_fine, alpha1M_av_times, label='$\\alpha_{\\mathit{1}avM}$', marker='o', markersize=2, color='tab:red')
ax5.set_ylabel('$\\alpha_{\\mathit{1}}$')
ax5.set_xlabel('$t, s$')
ax5.legend(loc=0)

# 7th set
ax6 = axs[2, 1]
ax6.plot(ts_out_fine, P0M_av_times, label='$P_{\\mathit{0}avM}$', marker='o', markersize=2, color='tab:cyan')
ax6.set_ylabel('$P_{\\mathit{0}}, Pa$')
ax6.set_xlabel('$t, s$')
ax6.legend(loc=0)

# 8th set
ax7 = axs[3, 1]
ax7.plot(ts_out_fine, uM_av_times, label='$u_{avM}$', marker='o', markersize=2, color='tab:orange')
ax7.set_ylabel('$u, m/s$')
ax7.set_xlabel('$t, s$')
ax7.legend(loc=0)

plt.tight_layout()
plt.show()
