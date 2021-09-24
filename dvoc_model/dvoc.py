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

        self.p = 0 
        self.q = 0

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
        # power at the terminal
        v = AlphaBeta.from_polar(self.v, self.theta)
        p = THREE_HALVES * (v.alpha * i.alpha + v.beta * i.beta)
        q = THREE_HALVES * (v.beta * i.alpha - v.alpha * i.beta)

        # feedback terms
        p_err = p - self.p
        q_err = q - self.q
        v_err = self.sin_phi * q_err + self.cos_phi * p_err
        theta_err = self.sin_phi * p_err - self.cos_phi * q_err


        dvdt = self.k0 * self.v * (self.two_v_nom_srd - 2*self.v**2) \
               - self.k1 * v_err / self.v

        dthetadt = self.omega_nom \
               - self.k1 * theta_err / (self.v**2)


        # %terminal voltage in alpha/beta frame
        # v_a = sqrt(2)*x(1)*cos(x(2));
        # v_b = sqrt(2)*x(1)*sin(x(2));

        # %power at the terminal
        # P = 3/2*(v_a*x(3) + v_b*x(4));
        # Q = 3/2*(v_b*x(3) - v_a*x(4));

        # %grid voltage
        # vg_a = sqrt(2)*vg_mag*cos(omega_nom*t);
        # vg_b = sqrt(2)*vg_mag*sin(omega_nom*t);

        # dx(1) = eps/k_v^2*x(1)*(2*V_nom^2 - 2*x(1)^2) - k_v*k_i/(3*C*x(1))*(sin(varphi)*(Q - Q_ref) + cos(varphi)*(P - P_ref)); %V
        # dx(2) = omega_nom - k_v*k_i/(3*C*x(1)^2)*(sin(varphi)*(P - P_ref) - cos(varphi)*(Q - Q_ref)); %theta

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
    dt = 1 / 10e3
    ts = np.arange(0, 500e-3, dt)

    # grid parameters
    grid = Dq0(SQRT_2*80, 0, 0)
    grid_omega = TWO_PI * 60

    # create a step function for dispatch (3A to 6A)
    ps = 250 * np.ones(len(ts))
    ps[len(ts)//2:] = 500

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
    for p, t in zip(ps, ts):
        # update dispatch
        dvoc.p = p

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