function [peak_hqm] = objective(unknowns, bicycle, speed)
% OBJECTIVE - Returns the maximum value of the HQM transfer function
% magnitude if the system is closed loop stable or an order of magnitude
% larger value based on the largest closed loop poles are if the system is
% closed loop unstable.
%
% Syntax: peak_hqm = objective(unknowns, bicycle, speed)
%
% Inputs:
%   unknowns - Vector of the optimization parameters (c trail, w wheelbase,
%              lam steer axis tilt, IFyy front wheel rotational inertia),
%              size 4x1.
%   bicycle - A char for the bicycle name, e.g. 'Benchmark', 'Pista', etc.
%   speed - Scalar value for the design speed.
%
% Outputs:
%   peak_hqm - Maximum scalar value of the HQM magnitude.

display(sprintf(repmat('=', 1, 79)))
display('Values passed into the objective function:')
display(sprintf('c: %1.5f ', unknowns(1)))
display(sprintf('w: %1.5f ', unknowns(2)))
display(sprintf('lam: %1.5f ', unknowns(3)))
display(sprintf('IFyy: %1.5f ', unknowns(4)))

freqs = linspace(0.01, 40, 200);

par = par_text_to_struct(['parameters/' bicycle 'Par.txt']);

% NOTE : Here to try for fun to see if different bicycle designs are needed
% for riding in low gravity.
% TODO : These should be moved to "optimal_bicycle.m".
%par.g = 1.6;  % acceleration due to gravity on the moon
%par.g = 3.71;  % acceleration due to gravity on mars

par.c = unknowns(1);  % trail
par.w = unknowns(2);  % wheelbase
par.lam = unknowns(3); % steer axis tilt

% Given an existing wheel with a rotational inertia, mass, and radius that
% approximately adhere to I=m*r^2, make sure that this relationship holds
% for the optimal front wheel inertias and that we only have to add mass to
% the existing wheel while keeping the radius constant or reduce the radius
% while keeping mass constant. These should be physically realizable barring
% that the radius doesn't get too tiny.
IFyy_opt = unknowns(4);

if IFyy_opt < par.IFyy
    % keep the mass the same but reduce the wheel radius
    par.rF = sqrt(IFyy_opt ./ par.mF);
elseif IFyy_opt > par.IFyy
    % keep radius of the wheel the same but add mass
    par.mF = IFyy_opt ./ par.rF.^2;
end

display(sprintf('rF: %1.5f ', par.rF))
display(sprintf('mF: %1.5f ', par.mF))

par.IFyy = IFyy_opt; % front wheel rotational inertia

[A, B, C, D] = whipple_pull_force_abcd(par, speed);

data = generate_data(bicycle, speed, ...
                     'simulate', false, ...
                     'forceTransfer', {}, ...
                     'fullSystem', false, ...
                     'stateSpace', {A, B, C, D}, ...
                     'display', 0);

%roll_loop = minreal(tf(data.closedLoops.Phi.num, data.closedLoops.Phi.den));
lateral_dev_loop = minreal(tf(data.closedLoops.Y.num, data.closedLoops.Y.den));

if ~isstable(lateral_dev_loop)
    % TODO : It may be best to return NAN here as per the CMAES statement
    % "An easy way to implement a hard non-linear constraint is to return
    % NaN. Then, this function evaluation is not counted and a newly sampled
    % point is tried immediately." Right now it returns a large number
    % relative to the target low HQM values.
    peak_hqm = max(10, 100 * max(real(pole(lateral_dev_loop))));
else
    num = data.handlingMetric.num;
    den = data.handlingMetric.den;
    [mag, ~, ~] = bode(tf(num, den), freqs);
    peak_hqm = max(mag);
end

display('Value of the objective function:')
peak_hqm
display(sprintf(repmat('=', 1, 79)))
