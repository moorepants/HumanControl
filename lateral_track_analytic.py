#!/usr/bin/env python

from sympy import Symbol, symbols, Matrix, zeros

# build the bicycle state space
b = {}
delta, deltaDot, phiDot, phi, psi, yP, yQ = symbols(
    'delta, deltaDot, phiDot, phi, psi, yP, yQ')
tPhi, tDelta, fB, tDeltaDot = symbols('tPhi, tDelta, fB, tDeltaDot')
b['x'] = Matrix([
    [phi],
    [delta],
    [phiDot],
    [deltaDot],
    [psi],
    [yP],
])
b['u'] = Matrix([
    [tPhi],
    [tDelta],
    [fB],
])
b['y'] = Matrix([
    [phi],
    [delta],
    [phiDot],
    [psi],
    [yQ],
])

numStates = len(b['x'])
numInputs = len(b['u'])
numOutputs = len(b['y'])

# Make the A and B matrices generic as they are just essentially submatrices in
# the closed loop system, but the C matrix plays a role in controller equations
# because we are tracking the lateral deviation of the front wheel contact. In
# this case, set the appropriate entries that are ones and zeros so that the
# final equatios will be in a reasonably understandable form.
b['A'] = zeros(numStates, numStates)
b['A'][0, 2] = 1
b['A'][1, 3] = 1
b['A'][2, 0] = Symbol('a_phidd_phi')
b['A'][2, 1] = Symbol('a_phidd_del')
b['A'][2, 2] = Symbol('a_phidd_phid')
b['A'][2, 3] = Symbol('a_phidd_deld')
b['A'][3, 0] = Symbol('a_deldd_phi')
b['A'][3, 1] = Symbol('a_deldd_del')
b['A'][3, 2] = Symbol('a_deldd_phid')
b['A'][3, 3] = Symbol('a_deldd_deld')
b['A'][4, 1] = Symbol('a_psid_del')
b['A'][4, 3] = Symbol('a_psid_deld')
b['A'][5, 4] = Symbol('a_ypd_psi')

b['B'] = zeros(numStates, numInputs)
b['B'][2, 1] = Symbol('b_phidd_Tdel')
b['B'][2, 2] = Symbol('b_phidd_F')
b['B'][3, 1] = Symbol('b_deldd_Tdel')
b['B'][3, 2] = Symbol('b_deldd_F')

b['C'] = zeros(numOutputs, numStates)
b['C'][0, 0] = 1  # phi
b['C'][1, 1] = 1  # delta
b['C'][2, 2] = 1  # phiDot
b['C'][3, 4] = 1  # psi
# yQ depends on the rear wheel location, the heading angle and the steer angle
b['C'][4, 1] = Symbol('c_yq_del')
b['C'][4, 4] = Symbol('c_yq_psi')
b['C'][4, 5] = 1

b['D'] = Matrix(numOutputs, numInputs, lambda i, j: 0.)

# Write the controller output as a function of the gains and the commanded
# lateral deviation.
kDelta, kPhiDot, kPhi, kPsi, kYQ = symbols('kDelta, kPhiDot, kPhi, kPsi, kYQ')
deltac, phiDotc, phic, psic, yc = symbols('deltac, phiDotc, phic, psic, yc')
psic = kYQ * (yc - yQ)
phic = kPsi * (psic - psi)
phiDotc = kPhi * (phic - phi)
deltac = kPhiDot * (phiDotc - phiDot)
Unm = kDelta * (deltac - delta)
Unm = Unm.expand()
controller = Matrix([
    Unm.coeff(phi),
    Unm.coeff(delta),
    Unm.coeff(phiDot),
    Unm.coeff(psi),
    Unm.coeff(yQ),
    Unm.coeff(yc),
]).T

# neuromuscular state space
omega = Symbol('omega')
zeta = Symbol('zeta')
n = {}
n['A'] = Matrix([[0, 1],
                 [-omega**2, -2*zeta*omega]])
n['B'] = Matrix([[0],
                 [omega**2]])
n['C'] = Matrix([1, 0]).T
n['D'] = Matrix([0])
n['x'] = Matrix([[tDelta],
                 [tDeltaDot]])
n['u'] = Matrix([Unm])
n['y'] = Matrix([tDelta])

# Compute the differential equations for the closed loop system.
xDot = b['A']*b['x'] + b['B']*b['u']
xDot = xDot.col_join(
    n['A']*n['x'] + n['B']*controller*(b['C']*b['x']).col_join(
        Matrix([yc])))
y = (b['C']*b['x']).col_join(n['C']*n['x'])

# Build the closed loop state space matrices.
systemStates = [phi, delta, phiDot, deltaDot, psi, yP, tDelta, tDeltaDot]
systemInputs = [fB, yc]
systemOutputs = [phi, delta, phiDot, psi, yQ, tDelta]


def mat_coeff(equations, variables, i, j):
    c = equations[i].expand().coeff(variables[j])
    if c is None:
        c = 0
    return c


A = Matrix(len(systemStates), len(systemStates),
           lambda i, j: mat_coeff(xDot, systemStates, i, j))

B = Matrix(len(systemStates), len(systemInputs),
           lambda i, j: mat_coeff(xDot, systemInputs, i, j))

C = Matrix(len(systemOutputs), len(systemStates),
           lambda i, j: mat_coeff(y, systemStates, i, j))
