function [peak_hqm] = objective(unknowns)
% OBJECTIVE - Returns the maximum value of the HQM transfer function
% magnitude if the system is closed loop stable or a large value based on
% how large the closed loop poles are if the system is closed loop unstable.
%
% Syntax: peak_hqm = objective(unknowns)
%
% Inputs:
%   unknowns - Vector of the optimization parameters (c trail, w wheelbase,
%              lam steer axis tilt, IFyy front wheel rotational inertia),
%              size 4x1.
% Outputs:
%   peak_hqm - Maximum scalar value of the HQM magnitude.

% NOTE : This value has to be manually set to be the target speed for the
% optimal design.
% TODO : Move setting this into "optimal_bicycle.m".
speed = 4.0;

display('==========================================')

display('Values passed into the objective function:')
unknowns

freqs = linspace(0.01, 40, 200);

% NOTE : This has to be manually set and match that in "optimal_bicycle.m".
par = par_text_to_struct('parameters/PistaPar.txt');

% NOTE : Here to try for fun to see if different bicycle designs are needed
% for riding in low gravity.
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
else IFyy_opt > par.IFyy
    % keep radius of the wheel the same but add mass
    par.mF = IFyy_opt / par.rF.^2;
end

par.IFyy = IFyy_opt; % front wheel rotational inertia

[A, B, C, D] = whipple_pull_force_abcd(par, speed);

data = generate_data('Benchmark', speed, ...
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
    % point is tried immediately."
    peak_hqm = max(10, 100 * max(real(pole(lateral_dev_loop))));
else
    num = data.handlingMetric.num;
    den = data.handlingMetric.den;
    [mag, ~, ~] = bode(tf(num, den), freqs);
    peak_hqm = max(mag);
end

display('Value of the objective function:')
peak_hqm
