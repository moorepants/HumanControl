function create_ieee_paper_plots(data, rollData)
% Creates all of the figures for the IEEE paper.
%
% Parameters
% ----------
% data : structure
%   A structure contating the data from generate_data.m for all of the bicycles
%   and speeds for the IEEE paper.
% rollData : structure
%   The data for a single bicycle at a single speed with roll torque as the
%   input.

global goldenRatio
% used for figure width to height ratio
goldenRatio = (1 + sqrt(5)) / 2;

% create a plot directory if one doesn't already exist
if exist('plots/', 'dir') ~= 7
    mkdir('plots/')
end

% Define some linestyles and colors for each of the six bicycles
linestyles = {'-', '-', '--', ...
              '--', '-.', '-.'};
colors = {'k', ...
          [0.5, 0.5, 0.5], ...
          'k', ...
          [0.5, 0.5, 0.5], ...
          'k', ...
          [0.5, 0.5, 0.5]};

loop_shape_example(data.Benchmark.Medium, 'Steer')
loop_shape_example(rollData, 'Roll')
plot_io_roll(rollData, 'Distance')
plot_io_roll(rollData, 'Time')
open_loop_all_bikes(data, linestyles, colors)
handling_all_bikes(data, linestyles, colors)
path_plots(data, linestyles, colors)
var = {'delta', 'phi', 'psi', 'Tdelta'};
io = {'output', 'output', 'output', 'input'};
typ = {'Distance', 'Time'};
for i = 1:length(var)
    for j = 1:length(typ)
        plot_io(var{i}, io{i}, typ{j}, data, linestyles, colors)
    end
end
phase_portraits(data.Benchmark.Medium)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loop_shape_example(bikeData, input)
% Creates the example loop shaping for the bicycle at medium speed.
%
% Parameters
% ----------
% bikeData : structure
%   Contains data for a single bicycle at a single speed.
% input : string
%   'Steer' or 'Roll' depending on what input was used to control the bicycle.

global goldenRatio

% closed loop bode plots
figure()
figWidth = 5.0;
figHeight = figWidth / goldenRatio;
set(gcf, ...
    'Color', [1, 1, 1], ...
    'PaperOrientation', 'portrait', ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0, 0, figWidth, figHeight], ...
    'PaperSize', [figWidth, figHeight])

freq = {0.1, 20.0};

hold all

closedLoops = bikeData.closedLoops;

if strcmp(input, 'Steer')
    linestyles = {'', '', '-.', '--', '-', '-.', '--', '-'};
    % the closed delta loop
    num = closedLoops.Delta.num;
    den = closedLoops.Delta.den;
    bodeplot(tf(num, den), freq);
    % a typical neuromuscular model
    num = 2722.5;
    den = [1, 13.96, 311.85, 2722.5];
    bodeplot(tf(num, den), freq);
    whichLines = 5:-1:3;
elseif strcmp(input, 'Roll')
    linestyles = {'', '', '-', '-'};
    whichLines = 4:-1:2;
else
    error('Bad input, use Steer or Roll')
end

% the closed phi dot loop
num = closedLoops.PhiDot.num;
den = closedLoops.PhiDot.den;
closedBode = bodeplot(tf(num, den), freq);

hold off

% clean it up
opts = getoptions(closedBode);
opts.Title.String = 'Closed Loop Bode Diagrams';
opts.YLim = {[-45, 20], [-360, 90]};
opts.PhaseMatching = 'on';
opts.PhaseMatchingValue = 0;
setoptions(closedBode, opts)

% find all the lines in the current figure
lines = findobj(gcf, 'type', 'line');
for i = 3:length(lines)
    set(lines(i), 'LineStyle', linestyles{i}, ...
                  'Color', 'k', ...
                  'LineWidth', 2.0)
end

plotAxes = findobj(gcf, 'type', 'axes');
% make the tick labels smaller
set(plotAxes(1), 'Fontsize', 8)
set(plotAxes(2), 'Fontsize', 8)
if strcmp(input, 'Steer')
    legWords = {'$\delta$', '$\dot{\phi}$','Neuromuscular model from [27]'};
elseif strcmp(input, 'Roll')
    legWords = {'$\dot{\phi}$'};
end
closeLeg = legend(lines(whichLines), ...
                  legWords, ...
                  'Location', 'Southwest', ...
                  'Interpreter', 'Latex', ...
                  'Fontsize', 8);

filename = ['benchmark' input 'Closed'];
pathToFile = ['plots' filesep filename];
print(gcf, '-deps2', '-loose', [pathToFile '.eps'])
fixPSlinestyle([pathToFile '.eps'])

% open loop plots
figure()
set(gcf, ...
    'Color', [1, 1, 1], ...
    'PaperOrientation', 'portrait', ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0, 0, figWidth, figHeight], ...
    'PaperSize', [figWidth, figHeight])

openLoops = bikeData.openLoops;

hold all

num = openLoops.Phi.num;
den = openLoops.Phi.den;
bodeplot(tf(num, den), freq);

num = openLoops.Psi.num;
den = openLoops.Psi.den;
bodeplot(tf(num, den), freq);

num = openLoops.Y.num;
den = openLoops.Y.den;
openBode = bodeplot(tf(num, den), freq);

hold off

% clean it up
opts = getoptions(openBode);
opts.Title.String = 'Open Loop Bode Diagrams';
opts.YLim = {[-80, 20], [-540, -80]};
opts.PhaseMatching = 'on';
opts.PhaseMatchingValue = 0;
opts.Grid = 'on';
setoptions(openBode, opts)

% find all the lines in the current figure
lines = findobj(gcf, 'type', 'line');
linestyles = {'', '', '-.', '--', '-', '-.', '--', '-'};
for i = 3:length(lines)
    set(lines(i), 'LineStyle', linestyles{i}, ...
                  'Color', 'k', ...
                  'LineWidth', 2.0)
end

plotAxes = findobj(gcf, 'type', 'axes');
% make the tick labels smaller
set(plotAxes(1), 'Fontsize', 8)
set(plotAxes(2), 'Fontsize', 8)
closeLeg = legend(lines(8:-1:6), ...
                  {'$\phi$', '$\psi$','$y$'}, ...
                  'Location', 'Southwest', ...
                  'Interpreter', 'Latex');

filename = ['benchmark' input 'Open.eps'];
pathToFile = ['plots' filesep filename];
print(gcf, '-deps2', '-loose', pathToFile)
fixPSlinestyle(pathToFile)

% handling qualities plot
num = bikeData.handlingMetric.num;
den = bikeData.handlingMetric.den;
w = linspace(0.01, 20, 200);
[mag, phase, freq] = bode(tf(num, den), w);
figure()
set(gcf, ...
    'Color', [1, 1, 1], ...
    'PaperOrientation', 'portrait', ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0, 0, figWidth, figHeight], ...
    'PaperSize', [figWidth, figHeight])

hold on

metricLine = plot(freq, mag(:)', 'k-', 'Linewidth', 2.0);
level1 = line([0, 20], [5, 5]);
level2 = line([0, 20], [8, 8]);
set(level1, 'Color', 'k', 'Linestyle', '--', 'Linewidth', 2.0)
set(level2, 'Color', 'k', 'Linestyle', '--', 'Linewidth', 2.0)

ylim([0, 10]);

ylabel('Handling Quality Metric')
xlabel('Frequency (rad/sec)')
text(3, 3, 'Level 1')
text(3, 6.5, 'Level 2')
text(3, 9, 'Level 3')
box on

filename = ['benchmark' input 'Handling.eps'];
pathToFile = ['plots' filesep filename];
print(gcf, '-deps2', '-loose', pathToFile)
fixPSlinestyle(pathToFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function open_loop_all_bikes(data, linestyles, colors)
% Creates open loop Bode plots of all the bikes.

global goldenRatio

bikes = fieldnames(data)

figure()
figWidth = 4.0;
set(gcf, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, figWidth, figWidth / goldenRatio], ...
    'PaperSize', [figWidth, figWidth / goldenRatio])

freq = {0.1, 20.0};

hold all
for i = 2:length(bikes)
    num = data.(bikes{i}).Medium.openLoops.Phi.num;
    den = data.(bikes{i}).Medium.openLoops.Phi.den;
    openBode = bodeplot(tf(num, den), freq);
end
hold off

% clean it up
opts = getoptions(openBode);
opts.Title.String = '$\phi$ Open Loop Bode Diagrams at 5 m/s';
opts.Title.Interpreter = 'Latex';
opts.YLim = {[-30, 10], [-540, -90]};
opts.PhaseMatching = 'on';
opts.PhaseMatchingValue = 0;
setoptions(openBode, opts)

% find all the lines in the current figure
plotAxes = findobj(gcf, 'type', 'axes');
magLines = findobj(plotAxes(2), 'type', 'line');
phaseLines = findobj(plotAxes(1), 'type', 'line');

for i = 2:length(magLines)
    set(magLines(i), ...
        'LineStyle', linestyles{i - 1}, ...
        'Color', colors{i - 1}, ...
        'LineWidth', 1.0)
    set(phaseLines(i), ...
        'LineStyle', linestyles{i - 1}, ...
        'Color', colors{i - 1}, ...
        'LineWidth', 1.0)
end

plotAxes = findobj(gcf, 'type', 'axes');
closeLeg = legend(magLines(2:7), ...
                  bikes(2:7), ...
                  'Location', 'Southwest');
filename = 'openBode.eps';
pathToFile = ['plots' filesep filename];
print(pathToFile, '-depsc')
fixPSlinestyle(pathToFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handling_all_bikes(data, linestyles, colors)
% Creates handling quality metric for all bikes.

global goldenRatio

bikes = fieldnames(data);
figure()

figWidth = 4.0;
set(gcf, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, figWidth, figWidth / goldenRatio], ...
    'PaperSize', [figWidth, figWidth / goldenRatio])

w = linspace(0.01, 20, 200);
speedNames = fieldnames(data.Browser);
fillColors = {[0.98, 0.98, 0.98], [0.93, 0.93, 0.93], [0.85, 0.85, 0.85]};
hold all
for j = 1:length(speedNames)
    magnitudes = zeros(length(w), length(bikes) - 1);
    for i = 2:length(bikes)
        num = data.(bikes{i}).(speedNames{j}).handlingMetric.num;
        den = data.(bikes{i}).(speedNames{j}).handlingMetric.den;
        [mag, phase, freq] = bode(tf(num, den), w);
        magnitudes(:, i - 1) = mag(:)';
    end
    maxMag = max(magnitudes, [], 2);
    area(freq, maxMag, ...
         'Facecolor', fillColors{j}, ...
         'Edgecolor', fillColors{j})
end

for j = 1:length(speedNames)
    metricLines = zeros(length(bikes) - 1, 1);
    for i = 2:length(bikes)
        num = data.(bikes{i}).(speedNames{j}).handlingMetric.num;
        den = data.(bikes{i}).(speedNames{j}).handlingMetric.den;
        [mag, phase, freq] = bode(tf(num, den), w);
        metricLines(i - 1) = plot(freq, mag(:)', ...
                                 'Color', colors{i - 1}, ...
                                 'Linestyle', linestyles{i - 1}, ...
                                 'Linewidth', 2.0);
    end
end
hold off

legend([{'2.5 m/s', '5.0 m/s', '7.5 m/s'}, bikes(2:end)'])

ylim([0, 20]);
level1 = line([0, 20], [5, 5]);
level2 = line([0, 20], [8, 8]);
set(level1, 'Color', 'k', 'Linestyle', '--', 'Linewidth', 1.0)
set(level2, 'Color', 'k', 'Linestyle', '--', 'Linewidth', 1.0)
ylabel('Handling Quality Metric')
xlabel('Frequency (rad/sec)')
text(3.5, 4.25, 'Level 1')
text(3, 6.5, 'Level 2')
text(3, 9, 'Level 3')
box on
filename = 'handling.eps';
pathToFile = ['plots' filesep filename];
print(pathToFile, '-depsc')
fixPSlinestyle(pathToFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path_plots(data, linestyles, colors)
% Creates a plot of the path tracking for all bikes at all speeds.

global goldenRatio

bikes = fieldnames(data);
speedNames = fieldnames(data.Browser);

figure()
figWidth = 5.0;
set(gcf, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, figWidth, figWidth / goldenRatio], ...
    'PaperSize', [figWidth, figWidth / goldenRatio])

hold all
for j = 1:length(speedNames)
    time = data.(bikes{2}).(speedNames{j}).time;
    path = data.(bikes{2}).(speedNames{j}).path;
    speed = data.(bikes{2}).(speedNames{j}).speed;
    plot(time * speed, path * j, 'k-')
    for i = 2:length(bikes)
        x = data.(bikes{i}).(speedNames{j}).outputs(:, 17);
        y = data.(bikes{i}).(speedNames{j}).outputs(:, 18);
        plot(x, y * j, 'Linestyle', linestyles{i - 1}, 'Color', colors{i - 1})
    end
end
hold off
xlim([30 190])
box on
legend(['Path', bikes(2:end)'])
xlabel('Distance (m)')
ylabel('Lateral Deviation (m)')
filename = 'paths.eps';
pathToFile = ['plots' filesep filename];
print(pathToFile, '-depsc')
fixPSlinestyle(pathToFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_io(variable, io, xAxis, data, linestyles, colors)
% Creates a plot of the time histories of a particular output or input variable
% for three different speeds with either time or distance on the x axis.
%
% Parameters
% ----------
% variable : string
%   The name of the variable you'd like to plot.
% io : string
%   'input' for input and 'output' for output.
% data : structure
%   Data for a set of bicycles, the first being the benchmark bicycle.
% xAxis : string
%   'Distance' or 'Time' on the x axis.
% linestyles : cell array
%   An array of linestyle types, one for each bicycle.
% colors : cell array
%   An array of colors, one for each bicycle.

global goldenRatio

if strcmp(io, 'input')
    names = {'Tphi',
             'Tdelta',
             'F'};
    prettyNames = {'$T_\phi$',
                   '$T_\delta$',
                   '$F$'};
    units = {'(Nm)',
             '(Nm)',
             '(N)'};
elseif strcmp(io, 'output')
    names = {'xP',
             'yP',
             'psi',
             'phi',
             'thetaP',
             'thetaR',
             'delta',
             'thetaF',
             'xPDot',
             'ypDot',
             'psiDot',
             'phiDot',
             'thetaPDot',
             'thetaRDot',
             'deltaDot',
             'thetaFDot',
             'xQ',
             'yQ'};
    prettyNames = {'$x_P$',
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
    units = {'(m)',
             '(m)',
             '(rad)',
             '(rad)',
             '(rad)',
             '(rad)',
             '(rad)',
             '(rad)',
             '(m/sec)',
             '(m/sec)',
             '(rad/sec)',
             '(rad/sec)',
             '(rad/sec)',
             '(rad/sec)',
             '(rad/sec)',
             '(rad/sec)',
             '(m)',
             '(m)'};
else
    error('Please choose i or o')
end

index = find(ismember(names, variable) == 1);

bikes = fieldnames(data);
speedNames = fieldnames(data.Browser);

figure()
figWidth = 5.0;
figHeight = figWidth / goldenRatio;
set(gcf, ...
    'Color', [1, 1, 1], ...
    'PaperOrientation', 'portrait', ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0, 0, figWidth, figHeight], ...
    'PaperSize', [figWidth, figHeight])

for j = 1:length(speedNames)
    subplot(3, 1, j)
    hold all
    for i = 2:length(bikes)
        oneSpeed = data.(bikes{i}).(speedNames{j});
        time = oneSpeed.time;
        speed = oneSpeed.speed;
        distance = time * speed;
        history = oneSpeed.([io 's'])(:, index);
        if strcmp(xAxis, 'Distance')
            plot(distance, history, ...
                 'Linestyle', linestyles{i - 1}, ...
                 'Color', colors{i - 1})
        elseif strcmp(xAxis, 'Time')
            plot(time, history, ...
                 'Linestyle', linestyles{i - 1}, ...
                 'Color', colors{i - 1})
        else
            error('Choose Time or Distance, no other')
        end
    end
    if strcmp(xAxis, 'Distance')
        xlabel('Distance (m)')
        xlim([30 190])
    else
        xlabel('Time (s)')
        xlim([30 / speed, 190 / speed])
    end
    first = [prettyNames{index} ' ' units{index}];
    second = sprintf(' at %1.1f m/s', speed);
    ylabel([first second], 'Interpreter', 'Latex')
    box on
    hold off
end
plotAxes = findobj(gcf, 'type', 'axes');
legend(plotAxes(3), bikes(2:end))
filename = [variable xAxis '.eps'];
print(['plots' filesep filename], '-depsc')
fixPSlinestyle(['plots' filesep filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_io_roll(rollData, xAxis)

global goldenRatio

% closed loop bode plots
figure()
figWidth = 5.0;
figHeight = figWidth / goldenRatio;
set(gcf, ...
    'Color', [1, 1, 1], ...
    'PaperOrientation', 'portrait', ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0, 0, figWidth, figHeight], ...
    'PaperSize', [figWidth, figHeight])

speed = rollData.speed;
time = rollData.time;
path = rollData.path;
frontWheel = rollData.outputs(:, 18);
rollAngle = rollData.outputs(:, 4);
steerAngle = rollData.outputs(:, 7);
rollTorque = rollData.inputs(:, 1);

% plot the path
subplot(2, 1, 1)
hold all
if strcmp(xAxis, 'Distance')
    plot(speed * time, path, 'k--', 'Linewidth', 1.0)
    plot(speed * time, frontWheel, 'k-', 'Linewidth', 1.0)
    xlabel('Distance (m)')
    xlim([30, 150])
elseif strcmp(xAxis, 'Time')
    plot(time, path, 'k--', 'Linewidth', 1.0)
    plot(time, frontWheel, 'k-', 'Linewidth', 1.0)
    xlabel('Time (s)')
    xlim([30 / speed, 150 / speed])
else
    error('Bad xAxis, choose Distance or Time')
end
hold off
box on
ylabel('Lateral Deviation (m)')
ylim([-0.2, 2.2])
legend({'Path', '$x_Q$'}, ...
       'Interpreter', 'Latex', ...
       'Fontsize', 8)

subplot(2, 1, 2)
hold all
if strcmp(xAxis, 'Distance')
    plot(speed * time, rollAngle, 'k-', 'Linewidth', 1.0)
    [ax, h1, h2] = plotyy(speed * time, steerAngle, speed * time, rollTorque);
    xlabel('Distance (m)')
    xlim(ax(1), [30, 150])
    xlim(ax(2), [30, 150])
elseif strcmp(xAxis, 'Time')
    plot(time, rollAngle, 'k-', 'Linewidth', 1.0)
    [ax, h1, h2] = plotyy(time, steerAngle, time, rollTorque);
    xlabel('Time (s)')
    xlim(ax(1), [30 / speed, 150 / speed])
    xlim(ax(2), [30 / speed, 150 / speed])
else
    error('Bad xAxis, choose Distance or Time')
end
hold off
box on

set(get(ax(1), 'Ylabel'), ...
    'String', 'Angle (rad)', ...
    'Color', 'k')
set(get(ax(2), 'Ylabel'), ...
    'String', 'Torque (Nm)', ...
    'Color', 'k')
set(ax, 'YColor', 'k', 'Fontsize', 8)
set(h1, 'Linestyle', '--', 'Color', 'k', 'Linewidth', 1.0)
set(h2, 'Linestyle', '-.', 'Color', 'k', 'Linewidth', 1.0)
legend({'$\phi$', '$\delta$', '$T_\phi$'}, ...
       'Interpreter', 'Latex', ...
       'Fontsize', 8, ...
       'Location', 'Southeast')

filename = ['roll' xAxis '.eps'];
print(['plots' filesep filename], '-depsc', '-loose')
fixPSlinestyle(['plots' filesep filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phase_portraits(bikeData)
% Creates four phase portrait plots.

global goldenRatio

figure()
figWidth = 5.0;
set(gcf, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, figWidth, figWidth / goldenRatio], ...
    'PaperSize', [figWidth, figWidth / goldenRatio])

gainChanges = [0.8, 1, 1, 1, 1;
               1, 1.2, 1, 1, 1;
               1, 1, 1.2, 1, 1;
               1, 1, 1.2, 0.8, 0.8];
loopNames = {'kDelta', 'kPhiDot', 'kPhi', 'kPhi'};
xy = [7, 15;
      4, 12;
      4, 12;
      4, 12];
xySource = {'outputs', 'outputsDot', 'outputs', 'outputs'};
xlabels = {'$\delta$ (rad)',
           '$\dot{\phi}$ (rad/s)',
           '$\phi$ (rad)',
           '$\phi$ (rad)'};
ylabels = {'$\dot{\delta}$ (rad/s)',
           '$\ddot{\phi}$ (rad/s$^2$)',
           '$\dot{\phi}$ (rad/s)',
           '$\dot{\phi}$ (rad/s)'};
legends = {'$k_\delta$ = ', '$k_{\dot{\phi}}$ = ',
           '$k_\phi$ = ', '$k_\phi$ = '};
floatSpec = {'%1.1f', '%1.3f', '%1.1f', '%1.1f'};

for i = 1:length(loopNames)
    a = {};
    for j = 1:length(gainChanges)
        a{j} = sprintf('%1.1f', gainChanges(i, j));
    end
    display(['Calculating gains as ' a{1} '*kDelta, ' a{2} ...
             '*kPhiDot, ' a{3} '*kPhi, ' a{4} '*kPsi, ' a{5} '*kY.'])
    twentyPercent = generate_data('Benchmark', 5.0, 'Steer', 0, ...
                                  gainChanges(i, :));
    subplot(2, 2, i)
    hold on

    x = bikeData.(xySource{i})(:, xy(i, 1));
    y = bikeData.(xySource{i})(:, xy(i, 2));
    plot(x, y, 'k-', 'Linewidth', 1.0)

    x = twentyPercent.(xySource{i})(:, xy(i, 1));
    y = twentyPercent.(xySource{i})(:, xy(i, 2));
    plot(x, y, 'k--', 'Linewidth', 1.0)

    hold off

    box on
    axis tight
    xlabel(xlabels{i}, 'Interpreter', 'Latex')
    ylabel(ylabels{i}, 'Interpreter', 'Latex')
    leg1 = sprintf(floatSpec{i}, bikeData.modelPar.(loopNames{i}));
    leg2 = sprintf(floatSpec{i}, twentyPercent.modelPar.(loopNames{i}));
    legend({[legends{i} leg1], [legends{i} leg2]} , 'Interpreter', 'Latex')
end

% save the plot
filename = 'phasePortraits.eps';
print(['plots' filesep filename], '-depsc')
fixPSlinestyle(['plots' filesep filename])
