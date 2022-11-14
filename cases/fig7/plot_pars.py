import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../../'))
sys.path.append(os.path.join(current_path, '../../../netgrid'))
sys.path.append(os.path.join(current_path, '../../../raplea'))

import raplea as rl

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
# plt.style.use('dracula')


def plot_pars(pars, S0):
    kR0 = []
    kR1 = []
    kR1Linear = []
    psi = []

    kR0D = []
    kR1D = []
    PCD = []

    PC = []
    f = []
    fD = []

    for i in S0:
        kR0.append(pars.kR0(i))
        kR1.append(pars.kR1(i))
        kR1Linear.append(pars.kR1Linear(i))
        psi.append(pars.psi(i, i))

        kR0D.append(pars.kR0D(i))
        kR1D.append(pars.kR1D(i))
        PCD.append(pars.PCD(i))

        PC.append(pars.PC(i))
        f.append(pars.f(i, i))
        fD.append(pars.fD(i, i))

    fig_width = 15
    y_scale = 0.8
    fig, axs = plt.subplots(3, 1, figsize=(fig_width / 3, fig_width * y_scale))

    # 1st set

    ax0 = axs[0]
    lns01 = ax0.plot(S0, kR0, label='$k_{r0}$', color='tab:blue')
    lns02 = ax0.plot(S0, kR1, label='$k_{r1}$', color='tab:red')
    lns03 = ax0.plot(S0, kR1Linear, label='$k_{r1L}$', color='tab:cyan')
    lns04 = ax0.plot(S0, f, label='$f$', color='tab:orange')
    ax0.set_ylabel('$k_{r}$, $f$')
    ax0.set_xlabel('$S_{0}$')
    ax02 = ax0.twinx()
    lns05 = ax02.plot(S0, psi, label='$\\psi$', color='tab:green')
    ax02.set_ylabel('$\\psi$')

    # added these three lines
    lns0 = lns01 + lns02 + lns03 + lns04 + lns05
    labs0 = [l.get_label() for l in lns0]
    ax0.legend(lns0, labs0, loc=0)

    # 2nd set
    ax1 = axs[1]
    lns11 = ax1.plot(S0, kR0D, label='$k_{r0}\'$', color='tab:blue')
    lns12 = ax1.plot(S0, kR1D, label='$k_{r1}\'$', color='tab:red')
    ax1.set_ylabel('$k_{r}\'$')
    ax1.set_xlabel('$S_{0}$')
    ax12 = ax1.twinx()

    lns13 = ax1.plot(S0, fD, label='$f\'$', color='tab:green')
    ax12.set_ylabel('$f\'$')

    # added these three lines
    lns1 = lns11 + lns12 + lns13
    labs1 = [l.get_label() for l in lns1]
    ax1.legend(lns1, labs1, loc=0)


    # 3nd set
    ax2 = axs[2]
    lns21 = ax2.plot(S0, PC, label='$P_{c}$', color='tab:green')
    ax2.set_xlabel('$S_{0}$')
    ax2.set_ylabel('$P_{c}, Pa$')

    ax22 = ax2.twinx()
    lns22 = ax22.plot(S0, PCD, label='$P_{c}\'$', color='tab:blue')
    ax22.set_ylabel('$P_{c}\', Pa$')

    # added these three lines
    lns2 = lns21 + lns22
    labs2 = [l.get_label() for l in lns2]
    ax2.legend(lns2, labs2, loc=0)

    plt.tight_layout()
    plt.show()


def plot_langmuir(pars, P0):
    fig_width = 4.8
    y_scale = 0.9

    fig, ax = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                           tight_layout=True)

    alpha = []
    alphaD = []
    for i in P0:
        alpha.append(pars.alpha1(i))
        alphaD.append(pars.alpha1D(i))

    lns1 = ax.plot(P0, alpha, label='$\\alpha$', color='tab:blue')
    ax.set_xlabel('$P_{0}, Pa$')
    ax.set_ylabel('$\\alpha$')

    ax2 = ax.twinx()
    lns2 = ax2.plot(P0, alphaD, label='$\\alpha\'$', color='tab:green')
    ax2.set_ylabel('$\\alpha\'$')

    # added these three lines
    lns = lns1 + lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    plt.show()


if __name__ == '__main__':

    alphaCoeffs = [-2.19813E-26, 3.73110E-19, 2.50417E-12, 9.35591E-06, 1.96098E-01]
    # alphaCoeffs = [0, 0, 0, 2E-06, 5]

    pars = rl.Pars(S0R=0.1, S1R=0.1, M0=0.465, M1=0.86, MC=1000., TC=100, m0=1.709,
                   m1=2., mC=2., mu0=2.E-5, mu1=5.E-4, alphaCoeffs=alphaCoeffs,
                   sorptionType='', is0Wet=False, dt=1.)

    S0 = np.linspace(0, 1, 300, dtype=float)
    plot_pars(pars, S0)
    press = np.linspace(5e5, 5e6, 300, dtype=float)
    plot_langmuir(pars, press)
