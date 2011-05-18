function generate_data(bike, speed, gain)
% Generates data files for the human operator control model.
%
% Parameters
% ----------
% bike : string
%   The shortname of the bicycle model to use.
% speed : float
%   The speed of the bicycle.
% gain : float
%   A general gain multiplier for the model.

warning off

% load the bicycle parameters
pathToParFile = ['parameters' filesep bike 'Par.txt'];
par = par_text_to_struct(pathToParFile);

% calculate the A, B, C, and D matrices
[A, B, C, D] = whipple_pull_force_abcd(par, speed)
initialConditions = [0, 0, 0, 0, 0, 0, 0, 0, ...
                     0.5, -speed / par.rR, 0];
outputs = {'xP',
           'yP',
           'psi',
           'phi',
           'thetaP',
           'thetaR',
           'delta',
           'thetaF',
           'xPdot',
           'yPdot',
           'psidot',
           'phidot',
           'thetaPdot',
           'thetaRdot',
           'deltadot',
           'thetaFdot',
           'xQ',
           'yQ'}

% human neuromuscular system
neuroNum = 900;
neuroDen = [1, 2 * .707 * 30, 900];

% preview time delay
timeDelay = 2.75

% load the gains
try
    pathToGainFile = ['gains' filesep bike 'Gains.txt'];
    [kDelta, kPhiDot, kPhi, kPsi, kY] = load_gains(pathToGainFile, speed)
catch
    display('No Gains found, all set to unity')
    kDelta = 0.0;
    kPhiDot = 0.0;
    kPhi = 0.0;
    kPsi = 0.0;
    kY = 0.0;
end

options = simset('SrcWorkspace', 'current');
sim('WhippleModel.mdl', [], options)

yQ = y(:, 2) + par.w .* y(:, 3) - par.c .* y(:, 7) .* cos(par.lambda)

figure(1)
plot(t, yQ, t, y(:, 18))

figure(2)
% plot the wheel contact points
subplot(5, 1, 1)
plot(t, y(:, 1), ...
     t, y(:, 2), ...
     t, y(:, 17), ...
     t, y(:, 18))
xlim([0, 5])
legend(outputs{[1, 2, 17, 18]})

subplot(5, 1, 2)
plot(t, y(:, 3), ...
     t, y(:, 4), ...
     t, y(:, 5), ...
     t, y(:, 6), ...
     t, y(:, 7), ...
     t, y(:, 8))
xlim([0, 5])
legend(outputs{[3, 4, 5, 6, 7, 8]})

subplot(5, 1, 3)
plot(t, y(:, 9), ...
     t, y(:, 10))
xlim([0, 5])
legend(outputs{[9, 10]})

subplot(5, 1, 4)
plot(t, y(:, 11), ...
     t, y(:, 12), ...
     t, y(:, 13), ...
     t, y(:, 15))
xlim([0, 5])
legend(outputs{[11, 12, 13, 15]})

subplot(5, 1, 5)
plot(t, y(:, 14), ...
     t, y(:, 16))
xlim([0, 5])
legend(outputs{[14, 16]})

function [kDelta, kPhiDot, kPhi, kPsi, kY] = load_gains(path, speed)
% Returns the gains from a file for a particular speed.
%
% Parameters
% ----------
% path : string
%   Path to a csv file with the gains for one bike at various speeds.
% speed : float
%   Speed associated with the desired gains.
%
% Returns
% -------
% kDelta : float
%   The gain for the steer angle loop.
% kPhiDot : float
%   The gain for the roll rate loop.
% kPhi : float
%   The gain for the roll loop.
% kPsi : float
%   The gain for the heading loop.
% kY : float
%   The gain for the lateral deviation loop.

contents = importdata(path);

% find the row with the corresponding speed
speeds = contents.data(:, 1);
row = find(speeds == speed);
if isempty(row)
    exception = MException('Eeeeeh:SpeedNotFound', ...
                           'No gain for speed %f', ...
                           speed);
    throw(exception)
else
    gains = contents.data(row, 2:end);
end

kDelta = gains(1);
kPhiDot = gains(2);
kPhi = gains(3);
kPsi = gains(4);
kY = gains(5);
