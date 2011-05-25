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
%
% Returns
% -------
% data : structure
%   Complete data output from the model. Includes:
%   - speed : float, speed of bicycle
%   - par : structure, bicycle parameters
%   - modelPar : structure, model input parameters
%   - closedLoops : structure, closed loop transfer functions
%   - openLoops : structure, open loop transfer functions
%   - time : vector (n, 1), time
%   - command : matrix (n, 5), commanded control
%   - inputs : matrix (n, 3), inputs to the bicycle system
%   - outputs : matrix (n, 18), outputs of the bicycle system

% this is for the latex expressions in the simulink model that can't compute
warning off

modelPar.gain = gain;

% load the bicycle parameters
pathToParFile = ['parameters' filesep bike 'Par.txt'];
par = par_text_to_struct(pathToParFile);
str = 'Parameters for the %s bicycle and rider have been loaded.';
display(sprintf(str, bike))

% calculate the A, B, C, and D matrices of the bicycle model
display(sprintf('Calculating the A, B, C, D matrices for %1.2f m/s', speed))
tic
[modelPar.A, modelPar.B, modelPar.C, modelPar.D] = ...
    whipple_pull_force_abcd(par, speed);
elapsedTime = toc;
display(sprintf('A, B, C, D calculated in %1.4f seconds.', elapsedTime))

% Keep in mind that the there is a function that relates steer angle, roll
% angle and pitch angle that must be enforced when setting any of those initial
% conditions.
modelPar.initialConditions = [-par.w, ... rear wheel contact x
                               0, ... rear wheel contact y
                               0, ... yaw angle
                               0, ... roll angle
                               0, ... pitch angle
                               0, ... rear wheel rotation
                               0, ... steer angle
                               0, ... front wheel rotation
                               0, ... roll rate
                               -speed / par.rR, ... rear wheel rate
                               0]; % steer rate

% human neuromuscular system
modelPar.neuroNum = 900;
modelPar.neuroDen = [1, 2 * 0.707 * 30, 900];

% handling qualities metric filter
modelPar.handlingFilterNum = 400;
modelPar.handlingFilterDen = [1, 40, 400];

% path filter
modelPar.pathFilterNum = (2.4 * gain)^2;
modelPar.pathFilterDen = [1, 2 * 2.4 * gain  (2.4*gain)^2];

% preview time delay
modelPar.timeDelay = 2.75;

% load the gains, set to zero if gains aren't available
try
    pathToGainFile = ['gains' filesep bike 'Gains.txt'];
    [modelPar.kDelta, modelPar.kPhiDot, modelPar.kPhi, ...
     modelPar.kPsi, modelPar.kY] = load_gains(pathToGainFile, speed);
catch
    display('No Gains found, all set to zero. No feedback.')
    modelPar.kDelta = 0.0;
    % the following two are just a hack to get around the 1/kPhi and 1/kPhiDot,
    % this can be handled better with some logic
    modelPar.kPhiDot = 1E-10;
    modelPar.kPhi = 1E-10;
    modelPar.kPsi = 0.0;
    modelPar.kY = 0.0;
end

% make a truth table for perturbing the loops
% the first row is default setup
perturbTable = [zeros(1, 5); eye(5)];
% make a truth table for closing the loops
closedTable = ~perturbTable;

% set all loops closed, perturbing none
modelPar.perturb = perturbTable(1, :);
modelPar.closed = closedTable(1, :);
modelPar.isHandling = 0;

loopNames = {'Delta', 'PhiDot', 'Phi', 'Psi', 'Y'};

% get the transfer functions for the closed loops
for i = 1:length(loopNames)
    str = 'Finding the closed loop transfer function of the %s loop.';
    display(sprintf(str, loopNames{i}))
    modelPar.loopNumber = i;
    modelPar.perturb = perturbTable(i + 1, :);
    update_model_variables(modelPar)
    [num, den] = linmod('WhippleModel');
    closedLoops.(loopNames{i}).num = num;
    closedLoops.(loopNames{i}).den = den;
end

% get the transfer functions for the open loops
for i = 1:length(loopNames)
    str = 'Finding the open loop transfer function of the %s loop.';
    display(sprintf(str, loopNames{i}));
    modelPar.loopNumber = i;
    % open the appropriate loop
    modelPar.perturb = perturbTable(i + 1, :);
    modelPar.closed = closedTable(i + 1, :);
    update_model_variables(modelPar);
    [num, den] = linmod('WhippleModel');
    openLoops.(loopNames{i}).num = num;
    openLoops.(loopNames{i}).den = den;
end

% get the handling quality metric
display('Finding the handling quality metric.')
modelPar.isHandling = 1;
modelPar.loopNumber = 3;
modelPar.closed = [0, 0, 1, 1, 1];
modelPar.perturb = [0, 0, 1, 0, 0];
update_model_variables(modelPar);
[num, den] = linmod('WhippleModel');
handlingMetric.num = num;
handlingMetric.den = den;

% close all the loops and simulate
modelPar.loopNumber = 0;
modelPar.isHandling = 0;
modelPar.closed = closedTable(1, :);
modelPar.perturb = perturbTable(1, :);
update_model_variables(modelPar)
display('Simulating the tracking task.')
tic;
sim('WhippleModel.mdl')
elapsedTime = toc;
display(sprintf('Simulation finished in %f1.3 seconds.', elapsedTime))

% set the initial point of the front wheel ahead of the rear wheel by the
% wheelbase length
y(:, 17) = y(:, 17) + par.w;

% write data for export
data.speed = speed;
data.par = par;
data.modelPar = modelPar;
data.closedLoops = closedLoops;
data.openLoops = openLoops;
data.handlingMetric = handlingMetric;
data.time = t;
data.command = command;
data.inputs = u;
data.outputs = y;
data.path = yc;

% plot
if basicPlots
    display('Making basic plots.')
    figure(1)
    % go through each loop and plot the bode plot for the closed loops
    hold all
    for i = 1:length(loopNames)
        num = closedLoops.(loopNames{i}).num;
        den = closedLoops.(loopNames{i}).den;
        bode(tf(num, den), {0.1, 20.0})
    end
    legend(loopNames)
    hold off

    figure(2)
    % go through each loop and plot the bode plot
    hold all
    for i = 1:length(loopNames)
        num = openLoops.(loopNames{i}).num;
        den = openLoops.(loopNames{i}).den;
        bode(tf(num, den), {0.1, 20.0})
    end
    legend(loopNames)
    hold off

    figure(3)
    num = handlingMetric.num;
    den = handlingMetric.den;
    wl = linspace(0.01, 20, 200);
    [mag, phase, freq] = bode(tf(num, den), wl);
    plot(wl, mag(:)')

    outputPlot = plot_outputs(t, y, yc);

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

function outputPlot = plot_outputs(t, y, yc)
% Returns a plot of the model outputs.
%
% Parameters
% ----------
% t : matrix, size(n, 1)
%   The time vector.
% y : matrix, size(n, 18)
%   The outputs of the bicycle system.
% yc : matrix, size(n, 1)
%   The path that was tracked.
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
plot(y(:, 17), yc, ...
     y(:, 1), y(:, 2), ...
     y(:, 17), y(:, 18))
legend({'Path', 'Rear Wheel', 'Front Wheel'})

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

function [kDelta, kPhiDot, kPhi, kPsi, kY] = load_gains(pathToGains, speed)
% Returns the gains from a file for a particular speed.
%
% Parameters
% ----------
% pathToGains : string
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

contents = importdata(pathToGains);

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
