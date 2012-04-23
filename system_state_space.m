function sys = system_state_space(stype, bicycle, gains, neuro, inputs, outputs)
% function sys = system_state_space(stype, bicycle, gains, neuro, inputs, outputs)
%
% Returns the closed loop system state space of the Hess bicycle-rider
% system for heading or lateral deviation tracking with bicycle steer torque
% as the controlled input.
%
% stype : string
%   This should specify 'lateral' or 'heading' for lateral deviation
%   tracking and heading tracking closed loop systems, respectively.
% bicycle : ss
%   A continuous state space system describing the bicycle dynamics. The
%   following parameters must be set.
%   A : matrix, m x m
%       The state matrix of the bicycle system. The minimum set of states
%       for the 'lateral' system type must be [yP, psi, phi, delta, phiDot,
%       deltaDot] and the minimum set of states for the 'heading' type must
%       be [psi, phi, delta, phiDot, deltaDot] to allow for the calculation
%       of the necessary output variables.
%   B : matrix, m X n
%       The input matrix of the bicycle system. The input must at least
%       include `tDelta`.
%   C : matrix, o x m
%       The output matrix of the bicycle system. The outputs must include at
%       least yQ, psi, phi, phiDot, and delta.
%   D : matrix, o X n
%       The feed forward matrix of the bicycle system. This is always zero.
%   StateName : cell array of chars, m x 1
%       The names of the states in order.
%   InputName : cell array of chars, n x 1
%       The names of the inputs in order.
%   OutputName : cell array of chars, o x 1
%       The names of the outputs in order.
% gains : matrix, 1 x 5 or 1 x 4
%   The feedback gains [kDelta, kPhiDot, kPhi, kPsi, kYq]. kYq will be
%   ignored if provided for the 'heading' system type.
% neuro : double, 1 x 2
%   The neuromuscular frequency and damping ratio, [wnm, zetanm].
% inputs : cell array of chars, m + 2 x q
%   The names of the desired q system external inputs in order. If the array
%   is empty, then the system B matrix is empty. 'yQc' is an option for the
%   lateral system and 'yPsic' for the heading system type for lateral
%   deviation and reference tracking inputs.
% outputs : cell array of chars, p x m + 2
%   The names of the desired system outputs in order. These can be a subset
%   of the outputs in the bicycle outputs in addition to `tDelta` and
%   `tDeltaDot`.
%
% Returns
% -------
% sys : ss
%   The closed loop 'lateral' or 'heading' state space system.
%   StateName : cell array of chars, m + 2 x 1
%       The names of the desired system states in order. This will always be
%       the states for the bicycle plus two additional states `tDelta` and
%       `tDeltaDot`.
%   InputName : cell array of chars, n x 1
%       The names of the desired system inputs in order. These will always
%       be the inputs of the bicycle minus `tDelta` and plus `yc`, the
%       commanded lateral deviation of the front wheel.

% build the state matrix
% append the two additional state names
states = [bicycle.StateName; 'tDelta'; 'tDeltaDot'];
% add room for the two new states, tDelta and tDeltaDot
[mA, nA] = size(bicycle.A);
A = zeros(mA + 2, nA + 2);
% put the bicycle equations in the upper left corner
A(1:mA, 1:nA) = bicycle.A;
% put the steer torque B column in the new A matrix
A(1:mA, index('tDelta', states)) = bicycle.B(:, index('tDelta', ...
    bicycle.InputName));

% put the neuro block state space in the bottom right corner
wnm = neuro(1); % the natural frequency of the neuromuscular mode
zetanm = neuro(2); % damping ratio of the neuromuscular mode
%zeta = 0.707; % damping ratio of the neuromuscular mode
A(mA + 1:end, nA + 1:end) = [0, 1; -wnm^2, -2 * zetanm * wnm];

% build a bicycle C matrix that is minimal for the feedback loop
if strcmp(stype, 'lateral')
    minBicycleOutputs = {'psi', 'phi', 'delta', 'phiDot', 'yQ'};
    minBicycleStates = {'yP', 'psi', 'phi', 'delta', 'phiDot', 'deltaDot'};
elseif strcmp(stype, 'heading')
    minBicycleOutputs = {'psi', 'phi', 'delta', 'phiDot'};
    minBicycleStates = {'psi', 'phi', 'delta', 'phiDot', 'deltaDot'};
end
% copy the bicycle C matrix
minC = bicycle.C;
% delete the rows and columns for the extra states and outputs
rows = find(~ismember(bicycle.OutputName, minBicycleOutputs));
cols = find(~ismember(bicycle.StateName, minBicycleStates));
minC(rows, :) = [];
minC(:, cols) = [];

% build the tDeltaDot equation
tDeltaDotRow = index('tDeltaDot', states);
kDelta = gains(1);
kPhiDot = gains(2);
kPhi = gains(3);
kPsi = gains(4);
for i = 1:length(minBicycleStates)
    coef = ...
    -kDelta * wnm^2 * ...
    (minC(1, i) * kPhi * kPhiDot * kPsi +...
     minC(2, i) * kPhi * kPhiDot + ...
     minC(3, i) + ...
     minC(4, i) * kPhiDot);
    if strcmp(stype, 'lateral')
        kYQ = gains(5);
        coef = coef + -kDelta * wnm^2 * minC(5, i) * kPhi * kPhiDot * kPsi * kYQ;
    end
    A(tDeltaDotRow, index(minBicycleStates{i}, states)) = coef;
end

% build the input matrix
% remove tDelta as an input, add yQc or yPsic as the last input to the
% system and add two rows for the two new state equations
[mB, ~] = size(bicycle.B);
B = zeros(mB + 2, length(inputs));
% for each input that is in bicycle.InputName that isn't steer torque, add
% the column into B
% if yQc or psic is inputs, add it's column
for i = 1:length(inputs)
    if strcmp(inputs{i}, 'yQc') && strcmp(stype, 'lateral')
        B(:, i) = [zeros(mB + 1, 1); wnm^2 * prod(gains)];
    elseif strcmp(inputs{i}, 'psic') && strcmp(stype, 'heading')
        B(:, i) = [zeros(mB + 1, 1); wnm^2 * prod(gains(1:4))];
    else
        B(1:mB, i) = bicycle.B(:, index(inputs{i}, bicycle.InputName));
    end
end

% build the output matrix
[~, nC] = size(bicycle.C);
C = zeros(length(outputs), nA + 2);
for i = 1:length(outputs)
    if strcmp(outputs{i}, 'tDelta')
        row = strcmp(outputs, 'tDelta');
        col = strcmp(states, 'tDelta');
        C(row, col) = 1;
    elseif strcmp(outputs{i}, 'tDeltaDot')
        row = strcmp(outputs, 'tDeltaDot');
        col = strcmp(states, 'tDeltaDot');
        C(row, col) = 1;
    else
        % grab the row in the bicycle output matrix and put it in the system
        % output matrix
        C(i, 1:nC) = bicycle.C(find(strcmp(bicycle.OutputName, outputs{i})), :);
    end
end

% build the feed forward matrix, this is always zero
D = zeros(length(outputs), length(inputs));

sys = ss(A, B, C, D);
sys.StateName = states;
sys.InputName = inputs;
sys.OutputName = outputs;
if strcmp(stype, 'lateral')
    sys.Name = 'Lateral Deviation Tracking';
elseif strcmp(stype, 'heading')
    sys.Name = 'Heading Tracking';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function is probably unecessary, as logical indexing can be use
% C(strcmp(value, vector), :) for example
function i = index(value, vector)
i = find(strcmp(value, vector));
