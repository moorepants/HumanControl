function sys = system_state_space(bicycle, gains, neuro, inputs, outputs)
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
%   states : cell array of chars, m x 1
%       The names of the states in order.
%   inputs : cell array of chars, n x 1
%       The names of the inputs in order.
%   outputs : cell array of chars, o x 1
%       The names of the outputs in order.
% gains : matrix, 5 x 1
%   The feedback gains [kDelta, kPhiDot, kPhi, kPsi, kY].
% neuro : double
%   The neuromuscular frequency.
% inputs : cell array of chars, m + 2 x q
%   The names of the desired system inputs in order. If the array is empty,
%   then the system B matrix is empty.
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

% build the state matrix
% append the two additional state names
states = [bicycle.states, 'tDelta', 'tDeltaDot'];
% add room for the two new states, tDelta and tDeltaDot
[mA, nA] = size(bicycle.A);
A = zeros(mA + 2, nA + 2);
% put the bicycle equations in the upper left corner
A(1:mA, 1:nA) = bicycle.A;
% put the steer torque B column in the new A matrix
A(1:mA, index('tDelta', states)) = bicycle.B(:, index('tDelta', bicycle.inputs));

% put the neuro block state space in the bottom right corner
zeta = 0.707; % damping ratio of the neuromuscular mode
A(mA + 1:end, nA + 1:end) = [0, 1; -neuro^2, -2 * zeta * neuro];

% build a bicycle C matrix that is minimal for the feedback loop
minBicycleOutputs = {'psi', 'phi', 'delta', 'phiDot', 'yQ'};
minBicycleStates = {'yP', 'psi', 'phi', 'delta', 'phiDot', 'deltaDot'};
% copy the bicycle C matrix
minC = bicycle.C;
% delete the rows and columns for the extra states and outputs
rows = find(~ismember(bicycle.outputs, minBicycleOutputs));
cols = find(~ismember(bicycle.states, minBicycleStates));
minC(rows, :) = [];
minC(:, cols) = [];

% build the tDeltaDot equation
tDeltaDotRow = index('tDeltaDot', states);
kDelta = gains(1);
kPhiDot = gains(2);
kPhi = gains(3);
kPsi = gains(4);
kYQ = gains(5);
for i = 1:length(minBicycleStates)
    A(tDeltaDotRow, index(minBicycleStates{i}, states)) = ...
    -kDelta * neuro^2 * ...
    (minC(1, i) * kPhi * kPhiDot * kPsi +...
     minC(2, i) * kPhi * kPhiDot + ...
     minC(3, i) + ...
     minC(4, i) * kPhiDot + ...
     minC(5, i) * kPhi * kPhiDot * kPsi * kYQ);
end

% build the input matrix
% remove tDelta as an input, add yc as the last input to the system and add
% two rows for the two new state equations
[mB, nB] = size(bicycle.B);
B = zeros(mB + 2, length(inputs));
% for each input that is in bicycle.inputs that isn't steer torque, add the
% column into B
% if yc is inputs, add it's column
for i = 1:length(inputs)
    if strcmp(inputs{i}, 'yc')
        B(:, i) = [zeros(mB + 1, 1); neuro^2 * prod(gains)];
    else
        B(1:mB, i) = bicycle.B(:, index(inputs{i}, bicycle.inputs));
    end
end

% the steer torque delta dot equation yc coefficient
% put the other bicycle inputs other than the steer torque into the new B
% matrix

% build the output matrix
[mC, nC] = size(bicycle.C);
C = zeros(length(outputs), nA + 2);
for i = 1:length(outputs)
    if strcmp(outputs{i}, 'tDelta')
        row = find(strcmp(outputs, 'tDelta'));
        col = find(strcmp(states, 'tDelta'));
        C(row, col) = 1;
    elseif strcmp(outputs{i}, 'tDeltaDot')
        row = find(strcmp(outputs, 'tDeltaDot'));
        col = find(strcmp(states, 'tDeltaDot'));
        C(row, col) = 1;
    else
        % grab the row in the bicycle output matrix and put it in the system
        % output matrix
        C(i, 1:nC) = bicycle.C(find(strcmp(bicycle.outputs, outputs{i})), :);
    end
end

% build the feed forward matrix, this is always zero
D = zeros(length(outputs), length(inputs));

sys.A = A;
sys.B = B;
sys.C = C;
sys.D = D;
sys.states = states;
sys.inputs = inputs;
sys.outputs = outputs;

function i = index(value, vector)
    i = find(strcmp(value, vector));
