function data = generate_data(bike, speed, varargin)
% function data = generate_data(bike, speed, varargin)
% Generates data files for the human operator control model.
%
% Parameters
% ----------
% bike : string
%   The name of the bicycle model to use. This corresponds to a file in the
%   ./parameters directory named <bike>Par.txt.
% speed : float
%   The speed of the bicycle.
% varargin : pairs of strings and values
%   input : string, optional
%       'Steer' or 'Roll'. 'Steer' is the default.
%   gains : matrix (5, 1), optional
%       If gains are given this will manually override the search for the
%       optimal gains. In order kDelta, kPhiDot, kPhi, kPsi, kY.
%   gainMuls : matrix (5, 1), optional
%       General gain multipliers. The gains are applied starting at the inner
%       loop going out. [1, 1, 1, 1, 1] is the default.
%   laneType : string, optional
%       'single' or 'double' for a single or double lane change maneuver.
%       'double' is the default.
%   neuroFreq : float, optional
%       The neuromuscular frequency. The default is 30 rad/sec.
%   plot : boolean, optional
%       If 1 basic plots will be shown, if 0 no plots will be shown. 0 is the
%       default.
%   simulate : boolean, optional
%       Default is true. If true the simulation will be available in the
%       output.
%   loopTransfer : boolean, optional
%       Default is true. If true the open and closed loop transfer functions
%       will be available in the output.
%   handlingQuality : boolean, optional
%       Default is true. If true the handling quality metric will be
%       available in the output.
%   forceTransfer : cell array of strings, optional
%       The default is {'Delta', 'PhiDot', 'Phi', 'Psi', 'Y', 'Tdelta'(or
%       'Tphi')}. The output will contain the transfer functions from
%       lateral force to steer angle, roll rate, roll angle, yaw angle,
%       lateral deviation and steer torque (or roll torque if input is
%       'Roll'). If the array is empty, then none of the transfer functions
%       are computed.  You can also provide a subset of the available
%       transfer functions.
%   stateSpace : cell array, optional
%       This cell array should contain the state space matrices {A, B, C, D}
%       for the whipple pull force bicycle model. Be sure that the `bike`
%       and `speed` match this state space model. If it isn't specified, the
%       state space calculation happens inside generate_data. This option
%       allows the potentially time intensive calculation of the state space
%       to happen outside of generate data. But be careful with it because
%       the arguments `bike` and `speed` become redundant.
%   fullSystem : boolean, optional
%       If true the state space matrices for the entire system with lateral
%       force as the only input are returned.
%   display : boolean, optional
%       If true the function will display information to screen as the
%       function runs else it will display nothing. The default is true.
%
% Returns
% -------
% data : structure
%   Complete data set from the model and simulations:
%   closedLoops : structure
%       Closed loop transfer functions for each loop.
%   command : matrix (n, 5)
%       The commanded input to each loop.
%   forceTF : structure
%       Transfer functions from pull force to various outputs.
%   gainMuls : matrix (5, 1)
%       Multipliers for each gain.
%   handlingMetric : structure
%       Transfer function for the handling quality metric.
%   inputs : matrix (n, 3)
%       Inputs to the bicycle model.
%   modelPar : structure
%       Simulink model input variables.
%       A : matrix (11, 11)
%           The state matrix. Refer to the documentation in
%           whipple_pull_force_abdc.m for details.
%       B : matrix (11, 3)
%           The input matrix.
%       C : matrix
%           The output matrix.
%       D : matrix
%           The feed forward matrix.
%       speed : float
%           The forward speed of the bicycle.
%       track : vector
%           The lateral coordinates of the desired track.
%       stoptime : float
%           The final time of the simulation.
%       initialConditions : matrix (11, 1)
%           The initial conditions for the simulation.
%       neuroNum : float
%           The numerator of the neuromuscular transfer function.
%       neuroDen : matrix (1, 3)
%           The coefficients to the denominator of the neuromuscular transfer
%           function.
%       pathFilterNum : float
%           The numerator of the path filter transfer function.
%       pathFileterDen : matrix (1, 3)
%           The coefficients to the denominator of the path filter transfer
%           function.
%       handlingFilterNum : float
%           The numerator of the handling quality metric filter transfer function.
%       handlingFileterDen : matrix (1, 3)
%           The coefficients to the denominator of the handling quality
%           metric filter transfer function.
%   openLoops : structure
%       Open loop transfer functions.
%   outputs : matrix (n, 18)
%       Outputs of the bicycle system.
%   outputsDots : matrix (n, 18)
%       Time derivatives of the outputs.
%   par : structure
%       Bicycle parameters.
%   path : matrix (n, 1)
%       Time delay adjusted path.
%   speed : float
%       Speed of bicycle.
%   time : matrix (n, 1)
%       Time.
%   system : structure
%       A : matrix (11, 11)
%           The state matrix.
%       B : matrix (11, 1)
%           The input matrix.
%       C : matrix (11, 1)
%           The output matrix.
%       D : matrix
%           The feed forward matrix.
%
% Examples
% --------
% % Generate the data set for the Benchmark bicycle at 5.0 m/s with roll as the
% % input.
% >>data = generate_data('Benchmark', 5.0, 'input', 'Roll');
%
% % Generate the data set for the Fisher bicycle at 7.5 m/s with steer input
% % and show the graphs.
% >>data = generate_data('Fisher', 7.5, 'input', 'Steer', 'plot', 1);
%
% % Generate the data set for the Browser bicycle at 2.5 m/s with steer as an
% % input and multiply the five gains by various values and show the graphs.
% >>data = generate_data('Browser', 2.5, 'plot', 1, 'gainMuls', [1.1, 1.1, 0.9, 1.0, 0.8])
%
% % Generate the data set for the Bianchi Pista bicycle at 7.5 m/s with steer as the
% % input and a single lane change as the manuever.
% >>data = generate_data('Pista', 7.5, 'laneType', 'single');

global CURRENT_DIRECTORY
% get the directory which this m-file is in
S = dbstack('-completenames');
[CURRENT_DIRECTORY, ~, ~] = fileparts(S(1).file);

% there are some unconnected ports in the simulink model that send out warnings
warning off

%% parse function arguments
% set the defaults for the optional arguments
defaults.input = 'Steer';
defaults.laneType = 'double';
defaults.gainMuls = [1, 1, 1, 1, 1];
defaults.gains = [];
defaults.neuroFreq = 30;
defaults.plot = 0;
defaults.simulate = 1;
defaults.loopTransfer = 1;
defaults.handlingQuality = 1;
defaults.forceTransfer = {'Delta', 'PhiDot', 'Phi', 'Psi', 'Y', 'Tdelta'};
defaults.stateSpace = {};
defaults.fullSystem = 1;
defaults.display = 1;

% load in user supplied settings
if size(varargin, 2) >= 1
    userSettings = varargin_to_structure(varargin);
else
    userSettings = struct();
end

% combine the defaults with the user settings
settings = overwrite_settings(defaults, userSettings);

global PRINT_TO_SCREEN
if settings.display
    PRINT_TO_SCREEN = 1;
else
    PRINT_TO_SCREEN = 0;
end

% show the bike and speed on the screen
display_if(sprintf(repmat('-', 1, 79)))
display_if(sprintf('%s at %1.2f m/s.', bike, speed))
display_if(sprintf(repmat('-', 1, 79)))

%% set model parameters
[modelPar, startLoop, par] = ...
    set_initial_model_parameters(bike, speed, settings);

% the name of the loops starting with the inner loop
loopNames = {'Delta', 'PhiDot', 'Phi', 'Psi', 'Y'};

%% set the gains
% if the user did not supply the gains, try to calculate them
if isempty(settings.gains)
    % give a warning that the program doesn't work well for low speeds
    if speed < 2.5
        display_if(sprintf(repmat('*', 1, 76)))
        display_if('Warning')
        display_if(sprintf(repmat('*', 1, 76)))
        s = ['The speed, %1.2f m/s, is less than 2.5 m/s. The PhiDot ', ...
             'loop often has a\nhard time converging. It is suggested ', ...
             'to supply the gains manually for these\nlower speeds.'];
        display_if(sprintf(s, speed))
        display_if(sprintf(repmat('*', 1, 76)))
    end
    % load the gain guesses
    pathToGainFile = [CURRENT_DIRECTORY filesep 'gains' filesep bike settings.input 'Gains.txt'];
    [modelPar.kDelta, modelPar.kPhiDot, modelPar.kPhi, ...
     modelPar.kPsi, modelPar.kY] = lookup_gains(pathToGainFile, speed);
    % now calculate exact feedback gains using successive loop closure
    modelPar = exact_gains(startLoop, loopNames, modelPar, settings);
else % set the gains as the user specified
    modelPar.kDelta = settings.gains(1);
    modelPar.kPhiDot = settings.gains(2);
    modelPar.kPhi = settings.gains(3);
    modelPar.kPsi = settings.gains(4);
    modelPar.kY = settings.gains(5);
end

% scale the gains
k = {'kDelta', 'kPhiDot', 'kPhi', 'kPsi', 'kY'};
kString = '';
for i = 1:length(k)
    modelPar.(k{i}) = settings.gainMuls(i) * modelPar.(k{i});
    kString = [kString sprintf('%s = %1.3f\n                  ', ...
               k{i}, modelPar.(k{i}))];
end
display_if(['Gains are set to: ', strtrim(kString)])

%% store transfer function data
if settings.loopTransfer
    perturbTable = [zeros(1, 5) % do not perturb any loop
                    eye(5)]; % perturb each loop individually
    % all the loops are closed at first
    closedTable = ones(6, 5);

    if strcmp(settings.input, 'Roll')
        % don't feed back delta
        closedTable(:, 1) = zeros(6, 1);
    end

    update_model_variables(modelPar)

    % store the transfer functions for the closed loops
    for i = startLoop:length(loopNames)
        str = 'Finding the %s closed loop transfer function.';
        display_if(sprintf(str, loopNames{i}))
        modelPar.loopNumber = i;
        modelPar.perturb = perturbTable(i + 1, :);
        modelPar.closed = closedTable(i + 1, :);
        update_model_variables(modelPar)
        [num, den] = linmod_switch('WhippleModel');
        closedLoops.(loopNames{i}).num = num;
        closedLoops.(loopNames{i}).den = den;
    end

    % make a truth table for closing the loops sequentially
    closedTable = ~perturbTable;

    % don't feed back delta if looking at roll control
    if strcmp(settings.input, 'Roll')
        closedTable(:, 1) = zeros(6, 1);
    end

    % get the transfer functions for the open loops
    for i = startLoop:length(loopNames)
        str = 'Finding the %s open loop transfer function.';
        display_if(sprintf(str, loopNames{i}));
        modelPar.loopNumber = i;
        modelPar.perturb = perturbTable(i + 1, :);
        modelPar.closed = closedTable(i + 1, :);
        update_model_variables(modelPar);
        [num, den] = linmod_switch('WhippleModel');
        openLoops.(loopNames{i}).num = num;
        openLoops.(loopNames{i}).den = den;
    end

    % store the loop transfer functions
    data.closedLoops = closedLoops;
    data.openLoops = openLoops;

    if settings.plot
        display_if('Generating loop transfer plots.')
        figure()
        % go through each loop and plot the bode plot for the closed loops
        hold all
        for i = startLoop:length(loopNames)
            num = closedLoops.(loopNames{i}).num;
            den = closedLoops.(loopNames{i}).den;
            bode(tf(num, den), {0.1, 100.0})
        end
        legend(loopNames(startLoop:end))
        title('Closed loop transfer functions')
        hold off

        figure()
        % go through each loop and plot the open loop bode plot
        hold all
        for i = startLoop:length(loopNames)
            num = openLoops.(loopNames{i}).num;
            den = openLoops.(loopNames{i}).den;
            bode(tf(num, den), {0.1, 100.0})
        end
        legend(loopNames(startLoop:end))
        title('Open loop transfer functions')
        hold off
    end
end

%% get the handling quality metric
if settings.handlingQuality
    display_if('Finding the handling quality metric.')
    % the handling qualities must be calculated with the phi loop at 2 rad/sec
    % crossover
    if strcmp(settings.input, 'Roll')
        origkPhi = modelPar.kPhi;
        % find the gain needed to move the current phi loop to a crossover of 2
        num = openLoops.Phi.num;
        den = openLoops.Phi.den;
        w = logspace(-1, 2, 1000);
        [mag,phase] = bode(tf(num,den), w);
        mag = mag(:)';
        phase = phase(:)';
        % get the magnitude at the desired crossover frequency
        MagCO = interp1(w, mag, 2.0);
        % calculate the gain needed to get the desired crossover frequency
        modelPar.kPhi = 1 / MagCO * origkPhi;
    end

    modelPar.isHandling = 1;
    modelPar.loopNumber = 3;
    modelPar.closed = [0, 0, 1, 1, 1];
    modelPar.perturb = [0, 0, 1, 0, 0];

    update_model_variables(modelPar);

    [num, den] = linmod_switch('WhippleModel');
    handlingMetric.num = num;
    handlingMetric.den = den;

    % change the gain
    if strcmp(settings.input, 'Roll')
        modelPar.kPhi = origkPhi;
    end

    % store the handling quality metric
    data.handlingMetric = handlingMetric;
    if settings.plot
        display_if('Generating handling quality plot.')
        figure()
        num = handlingMetric.num;
        den = handlingMetric.den;
        wl = linspace(0.01, 20, 100);
        [mag, phase, freq] = bode(tf(num, den), wl);
        plot(wl, mag(:)')
        title('Handling quality metric')
    end
end

%% find the transfer functions from pull force to various outputs
if ~isempty(settings.forceTransfer)
    ftf = settings.forceTransfer;
    modelPar.isHandling = 0;
    modelPar.isPullPerturb = 1;
    modelPar.perturb = [0, 0, 0, 0, 0];
    modelPar.closed = [1, 1, 1, 1, 1];

    % handle the fact that the user can specify both the control input and
    % either pull force to steer torque or roll torque
    if strcmp(settings.input, 'Steer')
        % replace 'Tphi'
        if any(ismember(ftf, 'Tphi'))
            display_if(['You have specified Steer as the input so Tphi ' ...
                     'will be replaced with Tdelta'])
            ftf{find(ismember(ftf, 'Tphi')==1)} = 'Tdelta';
        end
    elseif strcmp(settings.input, 'Roll')
        % don't feed back delta
        modelPar.closed = [0, 1, 1, 1, 1];
        % replace 'Tdelta'
        if any(ismember(ftf, 'Tdelta'))
            display_if(['You have specified Roll as the input so Tdelta ' ...
                     'will be replaced with Tphi'])
            ftf{find(ismember(ftf, 'Tdelta')==1)} = 'Tphi';
        end
    end

    % calculate each transfer function
    for i = 1:length(ftf)
        display_if(sprintf(['Calculating the pull force to %s transfer ' ...
                        'function.'], ftf{i}))
        if strcmp(ftf{i}, 'Tdelta') || strcmp(ftf{i}, 'Tphi')
            % Tdelta is connected to the 0 port in the multiswitch
            modelPar.loopNumber = 0;
        else
            modelPar.loopNumber = find(ismember(loopNames, ftf{i})==1);
        end

        update_model_variables(modelPar)

        [num, den] = linmod_switch('WhippleModel');
        data.forceTF.(ftf{i}).num = num;
        data.forceTF.(ftf{i}).den = den;
    end

    if settings.plot
        figure()
        hold all
        for i = 1:length(ftf)
            bode(tf(data.forceTF.(ftf{i}).num, ...
                    data.forceTF.(ftf{i}).den))
            legend(ftf)
        end
        hold off
    end

end

% get the full system state space
display_if('Calculating the full system state space with lateral input.')
modelPar.isHandling = 0;
modelPar.isPullPerturb = 1;
modelPar.isFullSystem = 1;
modelPar.perturb = [0, 0, 0, 0, 0];
modelPar.closed = [1, 1, 1, 1, 1];
modelPar.loopNumber = 0;

update_model_variables(modelPar)

[A, B, C, D] = linmod('WhippleModel');
% check to see if the final system is unstable
if any(real(eig(A)) > 0)
    display(sprintf(repmat('*', 1, 76)))
    display('Warning')
    display(sprintf(repmat('*', 1, 76)))
    s = ['The system is not stable with these gains. If the ', ...
             'simulation completes, the\ndata may be invalid. ', ...
             'Please give better gain guesses or supply the gains\n', ...
             'manually.'];
    display(sprintf(s))
    display(sprintf(repmat('*', 1, 76)))
    roots(G.den{:});
    display(settings.gains)
else % if not unstable
    % only save gains if not user supplied and the neuro frequency is
    % default
    if isempty(settings.gains) && (settings.neuroFreq - 30) < 1e-8
        % write the gains to file if the system is stable
        pathToGainFile = [CURRENT_DIRECTORY filesep 'gains' filesep bike ...
            settings.input 'Gains.txt'];
        newGains = [modelPar.kDelta, modelPar.kPhiDot, modelPar.kPhi, ...
        modelPar.kPsi, modelPar.kY];
        write_gains(pathToGainFile, speed, newGains)
        display_if(sprintf('Gains written to %s', pathToGainFile))
    else
        display_if('Gains were not saved to file.')
    end
end

if settings.fullSystem
    data.system.A = A;
    data.system.B = B;
    data.system.C = C;
    data.system.D = D;
    data.system.outputs = {'xP',
                           'yP',
                           'psi',
                           'phi',
                           'thetaB',
                           'thetaR',
                           'delta',
                           'thetaF',
                           'xPDot',
                           'yPDot',
                           'psiDot',
                           'phiDot',
                           'thetaBDot',
                           'thetaRDot',
                           'deltaDot',
                           'thetaFDot',
                           'xQ',
                           'yQ',
                           'tDelta'};
end

%% Simulate the system
if settings.simulate
    % close all the loops and simulate
    modelPar.loopNumber = 0;
    modelPar.isHandling = 0;
    modelPar.isPullPerturb = 0;
    modelPar.isFullSystem = 0;
    modelPar.perturb = [0, 0, 0, 0, 0];
    modelPar.closed = [1, 1, 1, 1, 1];

    update_model_variables(modelPar)

    display_if('Simulating the tracking task.')
    tic;
    sim('WhippleModel.mdl')
    elapsedTime = toc;
    display_if(sprintf('Simulation finished in %1.3f seconds.', elapsedTime))

    % set the initial point of the front wheel ahead of the rear wheel by the
    % wheelbase length
    y(:, 17) = y(:, 17) + par.w;

    % store simulation data
    data.time = t;
    data.command = command;
    data.inputs = u;
    data.outputs = y;
    data.outputsDot = yDot;
    data.path = yc;

    if settings.plot
        display_if('Generating simulation plots.')

        outputPlot = plot_outputs(t, y, yc);

        figure()
        plot(t, u)
        title('Inputs')
        xlabel('Time [s]')
        legend({'Roll Torque', 'Steer Torque', 'Pull Force'})
    end
end

% write data for export
data.speed = speed;
data.par = par;
data.modelPar = modelPar;
data.bicycle.states = {'xP', 'yP', 'psi', 'phi', 'thetaB', 'thetaR', 'delta', ...
             'thetaF', 'phiDot', 'thetaRDot', 'deltaDot'};
data.bicycle.outputs = {'xP', 'yP', 'psi', 'phi', 'thetaB', 'thetaR', 'delta', ...
             'thetaF', 'xPDot', 'yPDot', 'psiDot', 'phiDot', ...
             'thetaBDot', 'thetaRDot', 'deltaDot', 'thetaFDot', 'xQ', 'yQ'};
data.bicycle.inputs = {'tPhi', 'tDelta', 'fB'};

display_if(sprintf('Done. \n'))

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

outputs = {'$x_P$',
           '$y_P$',
           '$\psi$',
           '$\phi$',
           '$\theta_B$',
           '$\theta_R$',
           '$\delta$',
           '$\theta_F$',
           '$\dot{x}_P$',
           '$\dot{y}_P$',
           '$\dot{\psi}$',
           '$\dot{\phi}$',
           '$\dot{\theta}_B$',
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
    set(leg, 'interpreter', 'latex')
end

function k = find_closed_gain(loop, guess)
% Returns the gain required for a 10dB resonant peak in the closed loop
% transfer function.
%
% Parameters
% ----------
% loop : str
%   The name of the loop.
% guess : float
%   An initial guess for the gain.

k = fzero(@(x) delta_mag_closed(x, loop), guess);

function k = find_open_gain(loop, input)
% Returns the gain needed to set the crossover frequency at a desired value.
%
% Parameters
% ----------
% loop : string
%   The name of the loop.
% input : string
%   Whether this is 'Steer' or 'Roll' torque input
%
% Returns
% -------
% k : float
%   The gain needed for the desired crossover frequency.
%
% Notes
% -----
% This function assumes that all model parameters are set in the base workspace
% for the 'WhippleModel'.

% set the gain for this loop to unity
assignin('base', ['k' loop], 1)
% get the transfer function
[num, den] = linmod_switch('WhippleModel');
w = logspace(-1,2,1000);
[mag,phase] = bode(tf(num,den), w);
mag = mag(:)';
phase = phase(:)';

% set the desired open loop crossover frequency
if strcmp(loop, 'Phi')
    if strcmp(input, 'Steer')
        wBW = 2.0;
    elseif strcmp(input, 'Roll')
        wBW = 1.5;
    end
elseif strcmp(loop, 'Psi')
    if strcmp(input, 'Steer')
        wBW = 1.0;
    elseif strcmp(input, 'Roll')
        wBW = 0.75;
    end
elseif strcmp(loop, 'Y')
    if strcmp(input, 'Steer')
        wBW = 0.5;
    elseif strcmp(input, 'Roll')
        wBW = 0.375;
    end
end

% get the magnitude at the desired crossover frequency
MagCO = interp1(w, mag, wBW);
% calculate the gain needed to get the desired crossover frequency
k = 1 / MagCO;

function delta = delta_mag_closed(gain, loop)
% Returns the difference between a 10db overshoot and the overshoot calculated
% with the given gain.
%
% Parameters
% ----------
% gain : float
%   The gain for the closed loop.
% loop : str
%   The name of the loop.
%
% Returns
% -------
% delta : float
%   The difference in the resonant peak height and the desired height of 10db.

% set the gain for this loop
assignin('base', ['k' loop], gain)

% get the closed loop transfer function
[num, den] = linmod_switch('WhippleModel');
w = logspace(-2, 2, 1000);

% check for stability
%G = tf(num, den);
%if any(real(roots(G.den{:})) > 0)
    %display_if('Loop is not stable with this gain.')
    %roots(G.den{:})
%end

% uncomment to show the bode diagram at each step
%figure(25)
%bode(tf(num, den), w)

[mag, phase] = bode(tf(num, den), w);
% rewrite mag and phase so the dimension is reduced
mag = mag(:)';
phase = phase(:)';
% get the maximum magitude and index, this is the peak of the neurmuscular mode
[magmax, iMaxMag] = max(mag);
% find lower cutoff of 2 rad/sec
%lowi = min(find(w > 2));
% truncate the magnitude and frequency below resonance
magtrunc = mag(1:iMaxMag);
wtrunc = w(1:iMaxMag);
% differentiate the truncated magnitude
dmag = [0 diff(magtrunc)];
% differentiate it again
ddmag = [0 diff(dmag)];
% find the maximum which is just left of the far right zero crossing
[maxDD, indMaxDD] = max(ddmag);
% find the zero crossing in ddmag just left of its peak, this should occur
% at the local minima of dmag that corresponds to the flat point in mag, but
% sometimes it doesn't cross the zero line so just choose a point a certain
% distance from the peak (4 rad/sec to the left of the peak)
for i = indMaxDD:-1:2
    if ddmag(i - 1) < 0 && ddmag(i) > 0
        iMinDmag = i;
        %display_if(sprintf(['Found a zero crossing in ddmag ' ...
                        %'at the inflection point at %f.'], ...
                        %wtrunc(iMinDmag)))
        break
    else
        iMinDmag = 1;
    end
end

if iMinDmag == 1
    %display_if('Did not find a zero crossing in ddmag at the inflection point.')
    if strcmp(loop, 'Delta')
        [tmp, iMinDmag] = min(dmag);
        %display_if(sprintf('Setting the Delta low point to %f', ...
                        %wtrunc(iMinDmag)))
    elseif strcmp(loop, 'PhiDot')
        wMid = wtrunc(end) - 4;
        wMid = 3;
        [tmp, iMinDmag] = min(abs(wMid - wtrunc));
        %display_if(sprintf('Setting the PhiDot low point to %f', ...
                        %wtrunc(iMinDmag)))
    end
end

% uncomment to show the derivative of the magnitude at each step
%figure(20)
%subplot(3, 1, 1)
%plot(wtrunc, magtrunc, wtrunc(iMinDmag), magtrunc(iMinDmag), 'o')
%subplot(3, 1, 2)
%plot(wtrunc, dmag, wtrunc(iMinDmag), dmag(iMinDmag), 'o')
%subplot(3, 1, 3)
%plot(wtrunc, ddmag, wtrunc(iMinDmag), ddmag(iMinDmag), 'o')
%grid on
%pause

% get the magnitude at the flat part on the slope
magmin = magtrunc(iMinDmag);
% the ratio of magnitude of the resonant peak to the valley just left of the
% peak
rmag = magmax / magmin;
% the difference in the magnitude and the desired 10db peak
delta = rmag - sqrt(10);  % 10 dB corresponds to a magnitude ratio of sqrt(10)

function modelPar = exact_gains(startLoop, loopNames, modelPar, settings)
% function modelPar = exact_gains(startLoop, loopNames, modelPar, settings)
% Finds the exact values for the gains for each loop given the guesses
% provided in modelPar.
%
% Parameters
% ----------
% startLoop : integer
%   Either 1 or 2. If 1, all loops are set and if the Delta loop is skipped.
% loopNames : cell array of strings
%   The names of the loops starting with the Delta loop.
% modelPar : structure
%   The complete model parametes with the gains set as initial guesses.
% settings : structure
%   The combined default and user supplied optional settings.
%
% Returns
% -------
% modelPar : structure
%   The complete model parameters with the exact gains.

% make truth tables for perturbing and closing the loops
% the first row is default setup
perturbTable = [zeros(1, 5) % do not perturb any loop
                eye(5)]; % perturb each loop individually
closedTable = [1 1 1 1 1 % all loops closed
               1 0 0 0 0 % delta loop closed
               1 1 0 0 0 % delta, phidot loops closed
               1 1 0 0 0 % delta, phidot loops closed
               1 1 1 0 0 % delta, phidot, phi loops closed
               1 1 1 1 0]; % delta, phidot, phi, psi loops closed

if strcmp(settings.input, 'Roll')
    % don't feed back delta
    closedTable(:, 1) = zeros(6, 1);
end

for i = startLoop:length(loopNames)
    guess = modelPar.(['k' loopNames{i}]);
    str = ['Finding the loop transfer function ' ...
           'of the %s loop with a start guess of %1.4f.'];
    display_if(sprintf(str, loopNames{i}, guess))
    % set the logic for this loop calculation
    modelPar.loopNumber = i;
    modelPar.perturb = perturbTable(i + 1, :);
    modelPar.closed = closedTable(i + 1, :);
    update_model_variables(modelPar);
    if i == 1 || i == 2
        gain = find_closed_gain(loopNames{i}, guess);
    elseif i == 3 || i == 4 || i == 5
        gain = find_open_gain(loopNames{i}, settings.input);
    else
        error(sprintf('%s loop not found', loopNames{i}))
    end
    modelPar.(['k' loopNames{i}]) = gain;
    str = '%s loop gain set to %1.4f.';
    display_if(sprintf(str, loopNames{i}, gain))
end

function [modelPar, startLoop, par] = ...
    set_initial_model_parameters(bike, speed, settings)
% function [modelPar, startLoop, par] = set_initial_model_parameters(bike, speed, settings)
% Sets the model parameters based on the user input.
%
% Parameters
% ----------
% bike : string
%   The name of the bicycle model to use. This corresponds to a file in the
%   ./Parameters directory named <bike>Par.txt.
% speed : float
%   The speed of the bicycle.
% settings : structure
%   The combine default and user supplied optional settings.
%
% Returns
% -------
% modelPar : structure
%   The initial parameters for the simulink model.
% startLoop : integer
%   Either 1 or 2 depending if the input is steer or roll.
% par : structure
%   The physical parameters of the bicycle and rider.

global CURRENT_DIRECTORY

% set the speed
modelPar.speed = speed;

% generate the path to track
[pathX, pathY, pathT] = lane_change(35, 2, 0.2, 250, speed, ...
                                    500, settings.laneType, 60);
modelPar.track = [pathT, pathY];
modelPar.stopTime = pathT(end);

% load the bicycle parameters
pathToParFile = [CURRENT_DIRECTORY filesep 'parameters' filesep bike 'Par.txt'];
par = par_text_to_struct(pathToParFile);
str = 'Parameters for the %s bicycle and rider have been loaded.';
display_if(sprintf(str, bike))

% calculate the A, B, C, and D matrices of the bicycle model
if isempty(settings.stateSpace)
    display_if(sprintf('Calculating the A, B, C, D matrices for %1.2f m/s', speed))
    tic
    [modelPar.A, modelPar.B, modelPar.C, modelPar.D] = ...
        whipple_pull_force_abcd(par, speed);
    elapsedTime = toc;
    display_if(sprintf('A, B, C, D calculated in %1.4f seconds.', elapsedTime))
else
    display_if('A, B, C, D matrices already supplied')
    mats = {'A', 'B', 'C', 'D'};
    for i = 1:4
        modelPar.(mats{i}) = settings.stateSpace{i};
    end
end

% Keep in mind that the there is a function that relates steer angle, roll
% angle and pitch angle that must be enforced when setting any three of those
% initial conditions.
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
wnm = settings.neuroFreq;
display_if(sprintf('Neuromuscular frequency set to: %1.1f', wnm))
modelPar.neuroNum = wnm^2;
modelPar.neuroDen = [1, 2 * 0.707 * wnm, wnm^2];

% path filter
modelPar.pathFilterNum = 2.4^2;
modelPar.pathFilterDen = [1, 2 * 2.4, 2.4^2];

% handling qualities metric filter
modelPar.handlingFilterNum = 400;
modelPar.handlingFilterDen = [1, 40, 400];

% handling quality calculation flag
modelPar.isHandling = 0;

% set the pull perturb flag
modelPar.isPullPerturb = 0;

% set the full system flag to default
modelPar.isFullSystem = 0;

% set parameters for steer or roll inputs
if strcmp(settings.input, 'Steer')
    % start at the delta loop
    startLoop = 1;
    % use steer torque control
    modelPar.isRollInput = 0;
    % preview time delay
    modelPar.timeDelay = 2.75;
elseif strcmp(settings.input, 'Roll')
    % start at the phiDot loop
    startLoop = 2;
    % use the roll torque control
    modelPar.isRollInput = 1;
    % preview time delay
    modelPar.timeDelay = 3.5;
else
    error('Choose Steer or Roll as the input')
end

function [num, den] = linmod_switch(model)
% Returns the numerator and denominator of the current transfer function
% specified by the input and output blocks. This function is a hack and only
% exists because the full system state space switch has one vector input and
% one scalar input, which forces the output to always be a vector. A better
% solution for this may be to remove the switches related to the output and
% simply have one output that always outputs all the variables, instead of
% select ones.

[num, den] = linmod(model);

num = num(1, :);

function display_if(string)
% Prints the string to screen only if the print to screen global is true.

global PRINT_TO_SCREEN

if PRINT_TO_SCREEN
    disp(string)
end
