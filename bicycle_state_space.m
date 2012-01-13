function bicycle = bicycle_state_space(bicycle, speed)
% function bicycle = bicycle_state_space(bicycle, speed)
%
% Returns the state space system of the Whipple model linearized about the
% nominal configuration and the supplied speed.
%
% Parameters
% ----------
% bicycle : char
%   The name of a bicycle in the parameters directory.
% speed : double
%   The forward speed at which to linearize the model about.
%
% Returns
% -------
% bicycle : ss
%   The state space model of the bicycle.

% get the directory which this m-file is in
S = dbstack('-completenames');
[CURRENT_DIRECTORY, ~, ~] = fileparts(S(1).file);

% load the paramaters
par = par_text_to_struct([CURRENT_DIRECTORY filesep 'parameters' ...
    filesep bicycle 'Par.txt']);

% generate the state space matrices
[A, B, C, D] = whipple_pull_force_abcd(par, speed);

% name the states, outputs and inputs
stateNames = {'xP',
              'yP',
              'psi',
              'phi',
              'thetaP',
              'thetaR',
              'delta',
              'thetaF',
              'phiDot',
              'thetaRDot',
              'deltaDot'};

outputNames = {'xP',
               'yP',
               'psi',
               'phi',
               'thetaP',
               'thetaR',
               'delta',
               'thetaF',
               'xPDot',
               'yPDot',
               'psiDot',
               'phiDot',
               'thetaPDot',
               'thetaRDot',
               'deltaDot',
               'thetaFDot',
               'xQ',
               'yQ'};

inputNames = {'tPhi',
              'tDelta',
              'fB'};

% build the ss structure
bicycle = ss(A, B, C, D, ...
             'StateName', stateNames, ...
             'OutputName', outputNames, ...
             'InputName', inputNames);
