import sys
import os
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rc
from alive_progress import alive_bar, config_handler

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../../'))
sys.path.append(os.path.join(current_path, '../../../netgrid'))
sys.path.append(os.path.join(current_path, '../../../raplea'))

# from plot_analyt import plot_analyt
# from plot_pars import plot_pars
import raplea as rl

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"


# plt.style.use('dracula')


def deltJump(pars, S, Sini):
    return abs(pars.fD(S, S) - (pars.f(S, S) - pars.f(Sini, Sini)) / (S - Sini))


def calc_analytical(pars, grid, times, u_in):
    # ==========================================================================================
    #                                       input                                              #
    # ==========================================================================================
    Sini = pars._S0R
    poro = grid._phi
    L = grid._L
    Xini = 0.
    Xmax = L

    # ==========================================================================================
    #                                       analytics                                          #
    # ==========================================================================================

    SFmax = 0
    Fcur = 0
    SF_array = np.linspace(pars._S0R, 1. - pars._S1R, 1000)
    for SFmaxCur in SF_array:
        if pars.fD(SFmaxCur, SFmaxCur) < Fcur:
            SFmax = SFmaxCur
            break
        else:
            Fcur = pars.fD(SFmaxCur, SFmaxCur)

    deltJmpCur = deltJump(pars, SFmax, Sini)
    Sjmp = SFmax
    Sjmp_array = np.linspace(SFmax, 1. - pars._S1R, 1000)
    for SjmpCurr in Sjmp_array:
        if SjmpCurr > pars._S0R and SjmpCurr < 1. - pars._S1R:
            if deltJmpCur > deltJump(pars, SjmpCurr, Sini):
                Sjmp = SjmpCurr
                deltJmpCur = deltJump(pars, SjmpCurr, Sini)

    seriesS = [[] for _ in range(len(times))]
    seriesX = [[] for _ in range(len(times))]
    Xjmp = []

    for i in range(len(times)):
        Xjmp.append(times[i] * u_in / poro * pars.fD(Sjmp, Sjmp) + Xini)
        for Scur in np.linspace(1. - pars._S1R, Sjmp, 1000):
            Xcur = times[i] * u_in / poro * pars.fD(Scur, Scur) + Xini
            if Xcur <= Xmax:
                seriesS[i].append(Scur)
                seriesX[i].append(Xcur)
        if Xjmp[i] < Xmax:
            X_array = np.linspace(Xjmp[i], Xmax, 100)
            for Xcur in X_array:
                seriesS[i].append(Sini)
                seriesX[i].append(Xcur)

    return seriesX, seriesS, times


def calc_numerical(pars, grid, times, u_in):
    N = grid._N
    L = grid._L
    A = 0.0005
    S0R = pars._S0R
    S1R = pars._S1R

    S0_ini = S0R
    P0_ini = 100

    S0_left = 1.
    P0_left = 500.

    S0_right = 0.
    P0_right = P0_ini

    seriesSNum = []

    dt_ini = 1.0 / 4.
    t_max = times[-1]
    dt_out_coarse = t_max / 5.
    dt_out_fine = t_max / 100
    isCEndEffect = False

    alpha1_ini = pars.alpha1(P0_ini)

    u = [np.zeros(N + 1, dtype=np.float64), np.zeros(N + 1, dtype=np.float64)]
    P0 = [np.full(N, P0_ini, dtype=np.float64), np.full(N, P0_ini, dtype=np.float64)]
    S0 = [np.full(N, S0_ini, dtype=np.float64), np.full(N, S0_ini, dtype=np.float64)]
    alpha1 = [np.full(N, alpha1_ini, dtype=np.float64), np.full(N, alpha1_ini, dtype=np.float64)]

    locP = rl.LocP(grid, pars, u, P0, S0, alpha1)
    conP = rl.ConP(grid, pars, u, P0, S0, alpha1, isCEndEffect)
    equationP = rl.Equation(locP, conP, P0, isCEndEffect)
    equationP._solver = "sparseLU"
    # equationP._solver = "biCGSTAB"
    # equationP._tolerance = 1.e-16
    equationP._bcTypes = {'left': 'Neumann', 'right': 'Dirichlet'}
    equationP._bc = {'left': u_in * A, 'right': P0_right}

    locS = rl.LocS(grid, pars, u, P0, S0, alpha1)
    conS = rl.ConS(grid, pars, u, P0, S0, alpha1, isCEndEffect)
    equationS = rl.Equation(locS, conS, S0, isCEndEffect)
    equationS._solver = "sparseLU"
    # equationP._solver = "biCGSTAB"
    # equationS._tolerance = 1.e-16
    equationS._bcTypes = {'left': 'Neumann', 'right': 'Neumann'}
    equationS._bc = {'left': S0_left, 'right': S0_right}

    def cfd_sorption(equationP, equationS):
        equationP.calculateCoefficients()
        equationP.fillMatrix()
        equationP.solve()
        equationP.calculateU()
        equationP.calculateAlpha1()
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

    ts_out_coarse = [0]
    ts_out_fine = [0]

    is_t_out_coarse = False
    is_t_out_fine = False

    is_t_last = False
    t_cur = 0.
    dt_cur = 0.
    t_out_cur_coarse = dt_out_coarse
    t_out_cur_fine = dt_out_fine
    while True:

        dt_cur = dt_ini
        is_t_out_coarse = False
        is_t_out_fine = False

        if t_cur + dt_cur >= t_out_cur_fine:
            dt_cur = t_out_cur_fine - t_cur
            t_out_cur_fine += dt_out_fine
            is_t_out_fine = True

        if t_cur + dt_cur == t_out_cur_coarse:
            t_out_cur_coarse += dt_out_coarse
            is_t_out_coarse = True
        elif t_cur + dt_cur > t_out_cur_coarse:
            dt_cur = t_out_cur_coarse - t_cur
            t_out_cur_coarse += dt_out_coarse
            is_t_out_coarse = True
            t_out_cur_fine -= dt_out_fine
            is_t_out_fine = False

        if t_cur + dt_cur >= t_max:
            is_t_last = True
            is_t_out_coarse = True
            is_t_out_fine = True
            dt_cur = t_max - t_cur

        t_cur += dt_cur
        pars._dt = dt_cur
        cfd_sorption(equationP, equationS)

        iPrev = equationP._iPrev

        if is_t_out_coarse:
            ts_out_coarse.append(t_cur)
            S0_array_times.append(np.copy(S0[iPrev]))
            alpha1_array_times.append(np.copy(alpha1[iPrev]))
            P0_array_times.append(np.copy(P0[iPrev]))
            u_array_times.append(np.copy(u[iPrev]))

        if is_t_out_fine:
            ts_out_fine.append(t_cur)
            S0_av_times.append(np.mean(S0[iPrev]))
            alpha1_av_times.append(np.mean(alpha1[iPrev]))
            P0_av_times.append(np.mean(P0[iPrev]))
            u_av_times.append(np.mean(u[iPrev]))

        if is_t_last:
            break

    X_cells = np.linspace(0, grid._L, grid._N)

    return S0_array_times, X_cells


def plot_pars(pars, S0):
    kR0 = []
    kR1 = []

    for i in S0:
        kR0.append(pars.kR0(i))
        kR1.append(pars.kR1(i))

    fig_width = 3.5
    y_scale = 0.9
    fig, ax0 = plt.subplots(figsize=(fig_width, fig_width * y_scale))

    # 1st set
    lns01 = ax0.plot(S0, kR0, label='$k_{r\\mathit{0}}$', color='tab:blue')
    lns02 = ax0.plot(S0, kR1, label='$k_{r\\mathit{1}}$', color='tab:red')
    ax0.set_ylabel('$k_{r}$')
    ax0.set_xlabel('$S_{\\mathit{0}}$')

    # added these three lines
    lns0 = lns01 + lns02
    labs0 = [l.get_label() for l in lns0]
    ax0.legend(lns0, labs0, loc=0)

    plt.tight_layout()
    plt.savefig('fig6b_rel_perms.pdf', format="pdf",
                bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    grid = rl.Grid(phi=0.2, k=1E-12, L=0.1, A=0.0005, N=4000)
    alphaCoeffs = [0, 0, 0, 0, 0]
    pars = rl.Pars(S0R=0.1, S1R=0.1, M0=0.465, M1=0.86, MC=0., m0=1.709,
                   m1=2., mC=1., mu0=2.E-5, mu1=5.E-4, alphaCoeffs=alphaCoeffs,
                   sorptionType='', is0Wet=True, dt=0)

    S0 = np.linspace(0, 1, 300, dtype=float)
    plot_pars(pars, S0)

    time = 2000.
    u_in = 3.2E-6
    times = np.linspace(0, time, 6)
    seriesX, seriesS, times = calc_analytical(pars, grid, times, u_in)
    seriesSNum, X_cells = calc_numerical(pars, grid, times, u_in)

    # ==========================================================================================
    #                                       plotting                                          #
    # ==========================================================================================
    fig_width = 3.5
    y_scale = 0.9
    fig, ax = plt.subplots(figsize=(fig_width, fig_width * y_scale))
    ax.set_xlabel('$X, m$')
    ax.set_ylabel('$S_{\mathit{0}}$')

    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
              'tab:purple', 'tab:gray', 'tab:pink', 'tab:olive',
              'tab:cyan']

    line_names = ['analytical', 'numerical']
    line_types = ['-', '--']
    line1, = ax.plot([], [], linestyle=line_types[0], c='k', label=line_names[0])
    line2, = ax.plot([], [], linestyle=line_types[1], c='k', label=line_names[1])
    legend_1 = ax.legend(handles=[line1, line2], scatterpoints=1, frameon=True, labelspacing=1, loc=9)
    ax.add_artist(legend_1)

    lns = []
    for i in range(len(times)):
        line = ax.plot(seriesX[i], seriesS[i], color=colors[i], label='$\\mathit{' + str(round(times[i] / 100.)) + '}$')
        lns += line
        ax.plot(X_cells, seriesSNum[i], color='k', ls='--', linewidth=1.0)
    labs = [l.get_label() for l in lns]

    plt.ylim(0., 1.)
    plt.xlim(0., grid._L)
    ax.legend(lns, labs, title="$t, \\mathit{10^{\\mathit{2}}}s$", loc=1)
    plt.tight_layout()

    plt.savefig('fig6a_BL_validation.pdf', format="pdf",
                bbox_inches='tight')
    plt.show()
