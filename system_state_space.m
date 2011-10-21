function sys = system_state_space(bicycle, gains, neuro, outputs)
% Returns the full system state space of the bicycle rider controller with
% steer torque as the controlled input.
%
% bicycle : structure
%   A : matrix, m x m
%       The state matrix of the bicycle system. The minimum set of states
%       must be yP, psi, phi, delta, phiDot, deltaDot to allow for the
%       calculation of the necessary output variables.
%   B : matrix, m X n
%       The input matrix of the bicycle system. The input must at least
%       include `tDelta`.
%   C : matrix, o x m
%       The output matrix of the bicycle system. The outputs must include at
%       least yQ, psi, phi, phiDot, and delta.
%   D : matrix, o X n
%       The feed forward matrix of the bicycle system. This is always zero.
%   x : cell array of chars, m x 1
%       The names of the states in order.
%   u : cell array of chars, n x 1
%       The names of the inputs in order.
%   y : cell array of chars, o x 1
%       The names of the outputs in order.
% gains : matrix, 5 x 1
%   The feedback gains [kDelta, kPhiDot, kPhi, kPsi, kY].
% neuro : double
%   The neuromuscular frequency.
% ouputs : cell array of chars, p x m + 2
%   The names of the desired system outputs in order. These can be a subset
%   of the outputs in the bicycle outputs in addition to `tDelta` and
%   `tDeltaDot`.
%
% Returns
% -------
% sys : structure
%   states : cell array of chars, m + 2 x 1
%       The names of the desired system states in order. This will always be
%       the states for the bicycle plus two additional states `tDelta` and
%       `tDeltaDot`.
%   inputs : cell array of chars, n x 1
%       The names of the desired system inputs in order. These will always
%       be the inputs of the bicycle minus `tDelta` and plus `yc`, the
%       commanded lateral deviation of the front wheel.

% the outputs must have at least yQ, psi, phi, phiDot, and delta this also
% implies that the

% build the state matrix
% append the two additional state names
states = [bicycle.x, 'tDelta', 'tDeltaDot']
% add room for the two new states, tDelta and tDeltaDot
[mA, nA] = size(bicycle.A)
A = zeros(mA + 2, nA + 2)
% put the bicycle equations in the upper left corner
A(1:mA, 1:nA) = bicycle.A
% put the neuro block state space in the bottom right corner
zeta = 0.707 % damping ratio of the neuromuscular mode
A(mA + 1:end, nA + 1:end) = [0, 1;-neuro, -2 * zeta * neuro]
% put the steer torque B column in the new A matrix
A(1:mA, find(strcmp(states, 'tDelta'))) = bicycle.B(:, find(ismember(bicycle.u, 'tDelta')))
% get the row corresponding to yQ in bicycle.C
yQrow = find(strcmp('yQ', bicycle.y))
yPCol = find(strcmp('yP', bicycle.x))
psiCol = find(strcmp('psi', bicycle.x))
phiCol = find(strcmp('phi', bicycle.x))
deltaCol = find(strcmp('delta', bicycle.x))
phiDotCol = find(strcmp('phiDot', bicycle.x))
deltaDotCol = find(strcmp('deltaDot', bicycle.x))
% coefficients with respect to the bicycle states of the steer torque double dot equation
tDeltaDDot.yP =  -neuro^2 * prod(gains) * bicycle.C(yQrow, yPCol)
tDeltaDDot.psi = -neuro^2 * prod(gains(1:4)) * (1 + gains(1) * bicycle.C(yQrow, psiCol))
tDeltaDDot.phi = -neuro^2 * prod(gains(1:3)) * (1 + prod(gains(4:5)) * bicycle.C(yQrow, phiCol))
tDeltaDDot.delta = -neuro^2 * gains(1) * (1 + prod(gains(2:5)) * bicycle.C(yQrow, deltaCol))
tDeltaDDot.phiDot = -neuro^2 * prod(gains(1:2)) * (1 + prod(gains(3:5)) * bicycle.C(yQrow, phiDotCol))
tDeltaDDot.deltaDot = -neuro^2 * prod(gains) * bicycle.C(yQrow, deltaDotCol)
% add the steer torque double dot equation
tddcoefs = fieldnames(tDeltaDDot)
for i = 1:length(tddcoefs)
    % get the column of the state
    colNum = find(strcmp(states, tddcoefs{i}))
    % put the coefficients in the last row
    A(end, colNum) = tDeltaDDot.(tddcoefs{i})
end

% build the input matrix
% remove tDelta as an input, add yc as the last input to the system and add
% two rows for the two new state equations
[mB, nB] = size(bicycle.B)
B = zeros(mB + 2, nB)
% the steer torque delta dot equation yc coefficient
B(end, end) = neuro^2 * prod(gains)
% put the other bicycle inputs other than the steer torque into the new B
% matrix
inputs = {}
j = 1
for i = 1:length(bicycle.u)
    % skip the steer torque input because it is included in the new A matrix
    if ~strcmp(bicycle.u{i}, 'tDelta')
        B(1:mB, j) = bicycle.B(:, i)
        inputs{j} = bicycle.u{i}
        j = j + 1
    end
end
inputs = [inputs 'yc']

% build the output matrix
[mC, nC] = size(bicycle.C)
C = zeros(length(outputs), nA + 2)
for i = 1:length(outputs)
    if strcmp(outputs{i}, 'tDelta')
        row = find(strcmp(outputs, 'tDelta'))
        col = find(strcmp(states, 'tDelta'))
        C(row, col) = 1;
    elseif strcmp(outputs{i}, 'tDeltaDot')
        row = find(strcmp(outputs, 'tDeltaDot'))
        col = find(strcmp(states, 'tDeltaDot'))
        C(row, col) = 1;
    else
        % grab the row in the bicycle output matrix and put it in the system
        % output matrix
        C(i, 1:nC) = bicycle.C(find(strcmp(bicycle.y, outputs{i})), :)
    end
end

% build the feed forward matrix, this is always zero
D = zeros(length(outputs), length(inputs))

sys.A = A;
sys.B = B;
sys.C = C;
sys.D = D;
sys.states = states;
sys.inputs = inputs;
sys.outputs = outputs;
