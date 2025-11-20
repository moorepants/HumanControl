#!/usr/bin/env python

from sympy import Symbol, symbols, Matrix, zeros

# build the bicycle state space
b = {}
delta, deltadot, phidot, phi, psi, yP, yq = symbols(
    'delta, deltadot, phidot, phi, psi, yP, yq')
tphi, tdelta, fB, tdeltadot = symbols('tphi, tdelta, fB, tdeltadot')
b['x'] = Matrix([
    [phi],
    [delta],
    [phidot],
    [deltadot],
    [psi],
    [yP],
])
b['u'] = Matrix([
    [tphi],
    [tdelta],
    [fB],
])
b['y'] = Matrix([
    [phi],
    [delta],
    [phidot],
    [psi],
    [yq],
])

num_states = len(b['x'])
num_inputs = len(b['u'])
num_outputs = len(b['y'])

# Make the A and B matrices generic as they are just essentially submatrices in
# the closed loop system, but the C matrix plays a role in controller equations
# because we are tracking the lateral deviation of the front wheel contact. In
# this case, set the appropriate entries that are ones and zeros so that the
# final equatios will be in a reasonably understandable form.
b['A'] = zeros(num_states, num_states)
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

b['B'] = zeros(num_states, num_inputs)
b['B'][2, 1] = Symbol('b_phidd_Tdel')
b['B'][2, 2] = Symbol('b_phidd_F')
b['B'][3, 1] = Symbol('b_deldd_Tdel')
b['B'][3, 2] = Symbol('b_deldd_F')

b['C'] = zeros(num_outputs, num_states)
b['C'][0, 0] = 1  # phi
b['C'][1, 1] = 1  # delta
b['C'][2, 2] = 1  # phidot
b['C'][3, 4] = 1  # psi
# yq depends on the rear wheel location, the heading angle and the steer angle
b['C'][4, 1] = Symbol('c_yq_del')
b['C'][4, 4] = Symbol('c_yq_psi')
b['C'][4, 5] = 1

b['D'] = Matrix(num_outputs, num_inputs, lambda i, j: 0.)

# Write the controller output as a function of the gains and the commanded
# lateral deviation.
kdelta, kphidot, kphi, kpsi, kyq = symbols('kdelta, kphidot, kphi, kpsi, kyq')
deltac, phidotc, phic, psic, yc = symbols('deltac, phidotc, phic, psic, yc')
psic = kyq*(yc - yq)
phic = kpsi*(psic - psi)
phidotc = kphi*(phic - phi)
deltac = kphidot*(phidotc - phidot)
Unm = kdelta*(deltac - delta)
Unm = Unm.expand()
controller = Matrix([
    Unm.coeff(phi),
    Unm.coeff(delta),
    Unm.coeff(phidot),
    Unm.coeff(psi),
    Unm.coeff(yq),
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
n['x'] = Matrix([[tdelta],
                 [tdeltadot]])
n['u'] = Matrix([Unm])
n['y'] = Matrix([tdelta])

# Compute the differential equations for the closed loop system.
xdot = b['A']*b['x'] + b['B']*b['u']
xdot = xdot.col_join(
    n['A']*n['x'] + n['B']*controller*(b['C']*b['x']).col_join(
        Matrix([yc])))
y = (b['C']*b['x']).col_join(n['C']*n['x'])

# Build the closed loop state space matrices.
system_states = [phi, delta, phidot, deltadot, psi, yP, tdelta, tdeltadot]
system_inputs = [fB, yc]
system_outputs = [phi, delta, phidot, psi, yq, tdelta]


def mat_coeff(equations, variables, i, j):
    c = equations[i].expand().coeff(variables[j])
    if c is None:
        c = 0
    return c


A = Matrix(len(system_states), len(system_states),
           lambda i, j: mat_coeff(xdot, system_states, i, j))

B = Matrix(len(system_states), len(system_inputs),
           lambda i, j: mat_coeff(xdot, system_inputs, i, j))

C = Matrix(len(system_outputs), len(system_states),
           lambda i, j: mat_coeff(y, system_states, i, j))
