#!/usr/bin/env python

from sympy import Symbol, symbols, Matrix, zeros

# build the bicycle state space
b = {}
delta, deltaDot, phiDot, phi, psi, yP, yQ =\
    symbols('delta, deltaDot, phiDot, phi, psi, yP, yQ')
tPhi, tDelta, fB, tDeltaDot = symbols('tPhi, tDelta, fB, tDeltaDot')
b['x'] = Matrix([[yP], [psi], [phi], [delta], [phiDot], [deltaDot]])
b['u'] = Matrix([[tPhi], [tDelta], [fB]])
b['y'] = Matrix([[psi], [phi], [delta], [phiDot], [yQ]])

numStates = len(b['x'])
numInputs = len(b['u'])
numOutputs = len(b['y'])

# Make the A and B matrices generic as they are just essentially submatrices in
# the closed loop system, but the C matrix plays a role in controller equations
# because we are tracking the lateral deviation of the front wheel contact. In
# this case, set the appropriate entries that are ones and zeros so that the
# final equatios will be in a reasonably understandable form.
b['A'] = Matrix(numStates, numStates,
        lambda i, j: Symbol('a' + str(i + 1) + str(j + 1)))

b['B'] = Matrix(numStates, numInputs,
        lambda i, j: Symbol('b' + str(i + 1) + str(j + 1)))

#b['C'] = Matrix(numOutputs, numStates,
        #lambda i, j: Symbol('c' + str(i + 1) + str(j + 1)))

b['C'] = zeros((numOutputs, numStates))
b['C'][0, 1] = 1 # psi
b['C'][1, 2] = 1 # phi
b['C'][2, 3] = 1 # delta
b['C'][3, 4] = 1 # phiDot
# yQ depends on the rear wheel location, the heading angle and the steer angle
b['C'][4, 0] = 1
c52, c54 = symbols('c52 c54')
b['C'][4, 1] = c52
b['C'][4, 3] = c54

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
controller = Matrix([Unm.coeff(psi), Unm.coeff(phi), Unm.coeff(delta),
    Unm.coeff(phiDot), Unm.coeff(yQ), Unm.coeff(yc)]).T

# neuromuscular state space
omega = Symbol('w')
zeta = Symbol('zeta')
n = {}
n['A'] = Matrix([[0, 1],
                 [-omega**2, -2 * zeta * omega]])
n['B'] = Matrix([[0],
                 [omega**2]])
n['C'] = Matrix([1, 0]).T
n['D'] = Matrix([0])
n['x'] = Matrix([[tDelta],
                 [tDeltaDot]])
n['u'] = Matrix([Unm])
n['y'] = Matrix([tDelta])

# Compute the differential equations for the closed loop system.
xDot = b['A'] * b['x'] + b['B'] * b['u']
xDot = xDot.col_join(n['A'] * n['x'] + n['B'] * controller * (b['C'] *
        b['x']).col_join(Matrix([yc])))

# Build the closed loop state space matrices.
systemStates = [yP, psi, phi, delta, phiDot, deltaDot, tDelta, tDeltaDot]
systemInputs = [tPhi, fB, yc]
systemOutputs = [psi, phi, delta, phiDot, yQ, tDelta]

def mat_coeff(equations, variables, i, j):
    c = equations[i].expand().coeff(variables[j])
    if c is None:
        c = 0
    return c

A = Matrix(len(systemStates), len(systemStates), lambda i, j: mat_coeff(xDot,
    systemStates, i, j))

B = Matrix(len(systemStates), len(systemInputs), lambda i, j: mat_coeff(xDot,
    systemInputs, i, j))

y = (b['C'] * b['x']).col_join(n['C'] * n['x'])

C = Matrix(len(systemOutputs), len(systemStates), lambda i, j: mat_coeff(y,
    systemStates, i, j))
