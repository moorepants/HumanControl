function data = generate_data(bike, speed, gain, basicPlots)
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
% basicPlots : boolean
%   If 1 basic plots will be shown, if 0 no plots will be shown.

warning off

modelPar.gain = gain;

% load the bicycle parameters
pathToParFile = ['parameters' filesep bike 'Par.txt'];
par = par_text_to_struct(pathToParFile);
display(sprintf('Parameters for the %s bicycle and rider have been loaded.', bike))

% calculate the A, B, C, and D matrices of the bicycle model
display(sprintf('Calculating the A, B, C, D matrices for %1.2f m/s', speed))
tic
[modelPar.A, modelPar.B, modelPar.C, modelPar.D] = ...
    whipple_pull_force_abcd(par, speed);
elapsedTime = toc;
display(sprintf('A, B, C, D calculated in %1.4f seconds.', elapsedTime))

% more thought on the initial conditions of the front and rear wheels is needed
modelPar.initialConditions = [0, 0, 0, 0, 0, 0, 0, 0, ...
                              0.0, -speed / par.rR, 0];

% human neuromuscular system
modelPar.neuroNum = 900;
modelPar.neuroDen = [1, 2 * .707 * 30, 900];

% preview time delay
modelPar.timeDelay = 2.75;

% load the gains
try
    pathToGainFile = ['gains' filesep bike 'Gains.txt'];
    [modelPar.kDelta, modelPar.kPhiDot, modelPar.kPhi, ...
     modelPar.kPsi, modelPar.kY] = load_gains(pathToGainFile, speed);
catch
    display('No Gains found, all set to zero. No feedback.')
    modelPar.kDelta = 0.0;
    modelPar.kPhiDot = 0.0;
    modelPar.kPhi = 0.0;
    modelPar.kPsi = 0.0;
    modelPar.kY = 0.0;
end

% make a truth table for opening the loops
% zeros on the diagonal
closedTable = ones(5) - eye(5);
% add the default state of all closed loops
closedTable = [ones(1, 5); closedTable];

% set closed to the default, all loops closed
modelPar.closed = closedTable(1, :);

loopNames = {'Delta', 'PhiDot', 'Phi', 'Psi', 'Y'};

% get the transfer functions for the closed loops
for i = 1:length(loopNames)
    display(sprintf('Finding the closed loop transfer function of the %s loop.', loopNames{i}))
    modelPar.loopNumber = i;
    update_model_variables(modelPar);
    [num, den] = linmod('WhippleModel');
    closedLoops.(loopNames{i}) = [num; den];
end

% get the transfer functions for the open loops
for i = 1:length(loopNames)
    display(sprintf('Finding the open loop transfer function of the %s loop.', loopNames{i}));
    modelPar.loopNumber = i;
    % open the appropriate loop
    modelPar.closed = closedTable(i + 1, :);
    update_model_variables(modelPar);
    [num, den] = linmod('WhippleModel');
    openLoops.(loopNames{i}) = [num; den];
end

% close all the loops and simulate
modelPar.loopNumber = 0;
modelPar.closed = closedTable(1, :);
update_model_variables(modelPar)
display('Simulating the tracking task')
sim('WhippleModel.mdl')
display('Simulation finished')

% write data for export
data.speed = speed;
data.par = par;
data.modelPar = modelPar;
data.closedLoops = closedLoops;
data.openLoops = openLoops;
data.time = t;
data.command = command;
data.inputs = u;
data.outputs = y;

% plot
if basicPlots
    display('Making basic plots.')
    figure(1)
    % go through each loop and plot the bode plot for the closed loops
    hold all
    for i = 1:length(loopNames)
        nd = closedLoops.(loopNames{i});
        num = nd(1, :);
        den = nd(2, :);
        bode(tf(num, den), {0.1, 20.0})
    end
    legend(loopNames)
    hold off

    figure(2)
    % go through each loop and plot the bode plot
    hold all
    for i = 1:length(loopNames)
        nd = openLoops.(loopNames{i});
        num = nd(1, :);
        den = nd(2, :);
        bode(tf(num, den), {0.1, 20.0})
    end
    legend(loopNames)
    hold off

    outputPlot = plot_outputs(t, y);

    figure()
    plot(t, u)
    display('Plotting finished.')
end

function update_model_variables(modelPar)
% Puts all the variables needed for the simulink model in to the base
% workspace. This is a hack because linmod has no way to operate inside a
% function.
%
% Parameters
% ----------
% modelPar : structure
%   A structure that contains a field for each unknown variable in the simulink
%   model.

modelParNames = fieldnames(modelPar);
for i = 1:length(modelParNames)
    assignin('base', modelParNames{i}, modelPar.(modelParNames{i}))
end

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
           '$y_Q$'};

outputPlot = figure();
% plot the wheel contact points
subplot(6, 1, 1)
plot(y(:, 1), y(:, 2), ...
     y(:, 17), y(:, 18))
legend({'Rear Wheel', 'Front Wheel'})

plt.angles = [3, 4, 5, 7];
plt.wheelAngles = [6, 8];
plt.contactRates = [9, 10];
plt.rates = [11, 12, 13, 15];
plt.wheelRates = [14, 16];

pltFields = fieldnames(plt);
numPlots = length(pltFields);

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
