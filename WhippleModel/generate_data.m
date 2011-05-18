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
                     0.0, -speed / par.rR, 0];
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

outputPlot = plot_outputs(t, y)

figure(2)
plot(t, u)

function outputPlot = plot_outputs(t, y)
% Returns a plot of the model outputs.
%
% Parameters
% ----------
% t : matrix, size(n, 1)
%   The time vector.
% y : matrix, size(n, 18)
%
% Returns
% -------
% outputPlot : figure
%   Plot of the outputs versus time.

outputs = {'$x_P$',
           '$y_P$',
           '$\psi$',
           '$\phi$',
           '$\theta_P$',
           '$\theta_R$',
           '$\delta$',
           '$\theta_F$',
           '$\dot{x}_P$',
           '$\dot{y}_P$',
           '$\dot{\psi}$',
           '$\dot{\phi}$',
           '$\dot{\theta}_P$',
           '$\dot{\theta}_R$',
           '$\dot{\delta}$',
           '$\dot{\theta}_F$',
           '$x_Q$',
           '$y_Q$'}

outputPlot = figure()
% plot the wheel contact points
subplot(6, 1, 1)
plot(y(:, 1), y(:, 2), ...
     y(:, 17), y(:, 18))
legend({'Rear Wheel', 'Front Wheel'})

plt.angles = [3, 4, 5, 7]
plt.wheelAngles = [6, 8]
plt.contactRates = [9, 10]
plt.rates = [11, 12, 13, 15]
plt.wheelRates = [14, 16]

pltFields = fieldnames(plt)
numPlots = length(pltFields)

for i = 1:numPlots
    subplot(numPlots + 1, 1, i + 1)
    hold all
    numbers = plt.(pltFields{i});
    for j = 1:length(numbers)
        plot(t, y(:, numbers(j)))
    end
    hold off
    leg = legend(outputs{plt.(pltFields{i})});
    %set(leg, 'interpreter', 'latex')
end

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
