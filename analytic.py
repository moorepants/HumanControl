from sympy import Symbol, symbols, Matrix

# build the bicycle state space


b = {}
delta, deltaDot, phiDot, phi, psi, yP, yQ = symbols('delta, deltaDot, phiDot, phi, psi, yP, yQ')
tPhi, tDelta, fB, tDeltaDot = symbols('tPhi, tDelta, fB, tDeltaDot')
b['x'] = Matrix([[yP], [psi], [phi], [delta], [phiDot], ['deltaDot']])
b['u'] = Matrix([[tPhi], [tDelta], [fB]])
b['y'] = Matrix([[psi], [phi], [delta], [phiDot], [yQ]])

numStates = len(b['x'])
numInputs = len(b['u'])
numOutputs = len(b['y'])

b['A'] = Matrix(numStates, numStates,
        lambda i, j: Symbol('a' + str(i + 1) + str(j + 1)))

b['B'] = Matrix(numStates, numInputs,
        lambda i, j: Symbol('b' + str(i + 1) + str(j + 1)))

b['C'] = Matrix(numOutputs, numStates,
        lambda i, j: Symbol('c' + str(i + 1) + str(j + 1)))

b['D'] = Matrix(numOutputs, numInputs, lambda i, j: 0.)

# the controller block
kDelta, kPhiDot, kPhi, kPsi, kYQ = symbols('kDelta, kPhiDot, kPhi, kPsi, kYQ')
deltac, phiDotc, phic, psic, yc = symbols('deltac, phiDotc, phic, psic, yc')
Unm = kDelta * (deltac - delta)
delta_c = kPhiDot * (phiDotc - phiDot)
phiDot_c = kPhi * (phic - phi)
phi_c = kPsi * (psic - psi)
psi_c = kYQ * (yc - yQ)

Unm = Unm.subs({deltac: delta_c})
Unm = Unm.subs({phiDotc: phiDot_c})
Unm = Unm.subs({phic: phi_c})
Unm = Unm.subs({psic: psi_c}).expand()

controller = Matrix([Unm.coeff(psi), Unm.coeff(phi), Unm.coeff(delta),
    Unm.coeff(phiDot), Unm.coeff(yQ), Unm.coeff(yc)]).T

# neuromuscular state space
omega = Symbol('w')
zeta = Symbol('zeta')
n = {}
n['A'] = Matrix([[0, 1], [-omega, -2 * zeta * omega]])
n['B'] = Matrix([[0], [omega**2]])
n['C'] = Matrix([1, 0]).T
n['D'] = Matrix([0])
n['x'] = Matrix([[tDelta], [tDeltaDot]])
n['u'] = Matrix([Unm])
n['y'] = Matrix([tDelta])

# build the whole system
xDot = b['A'] * b['x'] + b['B'] * b['u']
xDot = xDot.col_join(n['A'] * n['x'] + n['B'] * controller * (b['C'] *
        b['x']).col_join(Matrix([yc])))

systemStates = [yP, psi, phi, delta, phiDot, deltaDot, tDelta, tDeltaDot]
systemInputs = [tPhi, fB, yc]
systemOutputs = [psi, phi, delta, phiDot, yQ, tDelta]

# build system state matrix
A = Matrix(len(systemStates), len(systemStates), lambda i, j:
        xDot[i].expand().coeff(systemStates[j]))

B = Matrix(len(systemStates), len(systemInputs), lambda i, j:
        xDot[i].expand().coeff(systemInputs[j]))

y = (b['C'] * b['x']).col_join(n['C'] * n['x'])

C = Matrix(len(systemOutputs), len(systemStates), lambda i, j:
        y[i].expand().coeff(systemStates[j]))
