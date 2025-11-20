import sympy as sm
import numpy as np
from scipy.signal import lsim, lti
import matplotlib.pyplot as plt

from lateral_track_analytic import (A, B, C, system_outputs)

args = list(sm.ordered(list(A.free_symbols | B.free_symbols | C.free_symbols)))

eval_sys = sm.lambdify(args, (A, B, C))

# q'' = Minv(-(g*K0 + v**2*K2)*q - v*C1*q' + T + H*F))
# Bicycle 1 from Hess, Moore, Hubbard 2012
M = np.array([[106.0, 1.55], [1.55, 0.25]])
K0 = np.array([[-93.2, -1.76], [-1.76, -0.68]])
K2 = np.array([[0.0, 77.3], [0.0, 1.58]])
C1 = np.array([[0.0, 29.9], [-0.45, 1.08]])
v, g = 5.0, 9.81
c, w, lam = 0.068, 1.12, 0.4
# H from Jason in Moore 2012 (should be Bicycle 1 values to be perfectly
# correct)
H = np.array([[0.943], [0.011]])
# Bicycle plant state space
Minv = np.linalg.inv(M)
# x = phi, del, phid, deld
Ab = np.block([
    [np.zeros_like(M), np.eye(2)],
    [-Minv@(g*K0 + v**2*K2), -Minv*v@C1]
])
# u = Tphi, Tdel, F
Bb = np.block([
    [np.zeros_like(M), np.zeros((2, 1))],
    [Minv, Minv@H]
])

vals = np.array([
    Ab[3, 1],  # a_deldd_del,
    Ab[3, 3],  # a_deldd_deld,
    Ab[3, 0],  # a_deldd_phi,
    Ab[3, 2],  # a_deldd_phid,
    Ab[2, 1],  # a_phidd_del,
    Ab[2, 3],  # a_phidd_deld,
    Ab[2, 0],  # a_phidd_phi,
    Ab[2, 2],  # a_phidd_phid,
    v/w*np.cos(lam),  # a_psid_del,
    c/w*np.cos(lam),  # a_psid_deld,
    v,  # a_ypd_psi,
    Bb[3, 2],  # b_deldd_F,
    Bb[3, 1],  # b_deldd_Tdel,
    Bb[2, 2],  # b_phidd_F,
    Bb[2, 1],  # b_phidd_Tdel,
    -c*np.cos(lam),  # c_yq_del,
    w,  # c_yq_psi,
    48.0,  # kDelta,
    9.03,  # kPhi,
    -0.08,  # kPhiDot,
    0.161,  # kPsi,
    0.097,  # kYQ,
    30.0,  # omega
    0.707,  # zeta
])

A_cl, B_cl, C_cl = eval_sys(*vals)

sys = lti(A_cl, B_cl, C_cl, np.zeros((C_cl.shape[0], B_cl.shape[1])))
t = np.linspace(0.0, 5.0, num=100)
u = np.zeros((len(t), B_cl.shape[1]))
u[len(t)//2:, 0] = 20.0  # F
u[:, 1] = 0.2  # yc

t, y, x = lsim(sys, u, t)

fig, axes = plt.subplots(nrows=y.shape[1])
for yi, ax, lab in zip(y.T, axes, system_outputs):
    ax.plot(t, yi, label=lab)
    ax.legend()
plt.show()
