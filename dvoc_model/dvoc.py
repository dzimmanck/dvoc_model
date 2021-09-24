from math import pi, sin, cos

from dvoc_model.reference_frames import SinCos, Abc, Dq0, AlphaBeta
from dvoc_model.constants import *

class Dvoc:
    def __init__(self,
                 eps: float = 15.,
                 k_v: float = 80.,
                 k_i: float = 0.2,
                 v_nom: float = 80.,
                 hz_nom: float = 60,
                 varphi: float = pi / 2,
                 l: float = 26.268e-6,
                 c: float = 0.2679
                 ):

        self.v_nom = v_nom
        self.hz_nom = hz_nom
        # self.omega_nom = 1 / sqrt(l*c)
        self.omega_nom = TWO_PI * hz_nom

        self.ia = 0 
        self.ir = 0

        # initialize state variables
        self.v = SQRT_2*80
        self.theta = 0

        # calculate oscillator constants
        self.k0 = eps/(k_v**2)
        self.k1 = k_v * k_i / (3 * c)
        self.two_v_nom_srd = 2*self.v_nom**2
        self.sin_phi = sin(varphi)
        self.cos_phi = cos(varphi)



    def step(self, dt, i):
        # currents
        sin_cos = SinCos.from_theta(self.theta)
        idq0 = i.to_dq0(sin_cos)
        ia = idq0.d
        ir = -idq0.q

        # feedback terms
        ia_err = ia - self.ia
        ir_err = ir - self.ir
        v_err = self.sin_phi * ir_err + self.cos_phi * ia_err
        theta_err = self.sin_phi * ia_err - self.cos_phi * ir_err

        dvdt = self.k0 * self.v * (self.two_v_nom_srd - 2*self.v**2) \
               - self.k1 * v_err

        dthetadt = self.omega_nom \
               - self.k1 * theta_err / self.v

        self.v += dvdt * dt
        self.theta += dthetadt * dt



if __name__ == "__main__":
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as mp

    # filter values
    Lf = 1.5e-3
    Rf = 0.4

    # simulation parameters
    dt = 1 / 100e3
    ts = np.arange(0, 500e-3, dt)

    # grid parameters
    grid = Dq0(SQRT_2*80, 0, 0)
    grid_omega = TWO_PI * 60

    # create a step function for dispatch (3A to 6A)
    ias = 3 * np.ones(len(ts))
    ias[len(ts)//2:] = 6

    # create an oscillator using the defaults
    dvoc = Dvoc()

    # dictionary for containing simulation results
    data = {'v_a': [],
            'v_b': [],
            'v_c': [],
            'i_a': [],
            'i_b': [],
            'i_c': []}

    # run simulation
    i = AlphaBeta(0, 0, 0)  # start out current at 0A
    for ia, t in zip(ias, ts):
        # update dispatch
        dvoc.ia = ia

        # update the grid voltage
        sin_cos = SinCos.from_theta(grid_omega * t)
        vg = grid.to_alpha_beta(sin_cos)

        # update the virtual oscillator
        dvoc.step(dt, i)
        v = AlphaBeta.from_polar(dvoc.v, dvoc.theta)

        # simulate the currents
        i.alpha += (dt/Lf*(v.alpha - vg.alpha - Rf*i.alpha))
        i.beta += (dt/Lf*(v.beta - vg.beta - Rf*i.beta))

        # update the data
        v_abc = v.to_abc()
        i_abc = i.to_abc()
        data['v_a'].append(v_abc.a)
        data['v_b'].append(v_abc.b)
        data['v_c'].append(v_abc.c)
        data['i_a'].append(i_abc.a)
        data['i_b'].append(i_abc.b)
        data['i_c'].append(i_abc.c)

    # plot the results
    data = pd.DataFrame(index=ts, data=data)
    ax = data.plot(y='i_a')
    data.plot(y='i_b', ax=ax)
    data.plot(y='i_c', ax=ax)
    mp.show()