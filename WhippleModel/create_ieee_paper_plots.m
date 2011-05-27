function create_ieee_paper_plots(data)
% Creates all of the figures for the IEEE paper.
%
% Parameters
% ----------
% data : structure
%   A structure contating the data from generate_data.m for all of the bicycles
%   and speeds for the IEEE paper.

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

%loop_shape_example(data.Benchmark)
%open_loop_all_bikes(data, linestyles, colors)
%handling_all_bikes(data, linestyles, colors)
%path_plots(data, linestyles, colors)
%plot_io('delta', 'output', data, linestyles, colors)
%plot_io('phi', 'output', data, linestyles, colors)
%plot_io('psi', 'output', data, linestyles, colors)
%plot_io('Tdelta', 'input', data, linestyles, colors)
phase_portraits(data.Benchmark.Medium)

function loop_shape_example(bikeData)
% Creates the example loop shaping for the bicycle at medium speed.

% used for figure width to height ratio
goldenRatio = (1 + sqrt(5)) / 2;

figure()
freq = {0.1, 20.0};

hold all
% the closed delta loop
closedLoops = bikeData.Medium.closedLoops;
num = closedLoops.Delta.num;
den = closedLoops.Delta.den;
deltaBode = bodeplot(tf(num, den), freq);

% the closed phi dot loop
num = closedLoops.PhiDot.num;
den = closedLoops.PhiDot.den;
closedBode = bodeplot(tf(num, den), freq);

% a typical neuromuscular model
num = 2722.5;
den = [1, 13.96, 311.85, 2722.5];
closedBode = bodeplot(tf(num, den), freq);

% clean it up
opts = getoptions(closedBode);
opts.Title.String = 'Closed Loop Bode Diagrams';
opts.YLim = {[-45, 20], [-360, 90]};
opts.PhaseMatching = 'on';
opts.PhaseMatchingValue = 0;
setoptions(closedBode, opts)
hold off

% find all the lines in the current figure
lines = findobj(gcf, 'type', 'line');
linestyles = {'', '', '-.', '--', '-', '-.', '--', '-'};
for i = 3:length(lines)
    set(lines(i), 'LineStyle', linestyles{i}, ...
                  'Color', 'k', ...
                  'LineWidth', 2.0)
end

plotAxes = findobj(gcf, 'type', 'axes');
closeLeg = legend(lines(8:-1:6), ...
                  {'$\delta$', '$\dot{\phi}$','Neuromuscular'}, ...
                  'Location', 'Southeast', ...
                  'Interpreter', 'Latex');

figWidth = 3.0;
set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])

filename = 'benchmarkClosed.eps';
pathToFile = ['plots' filesep filename];
saveas(gcf, pathToFile)
fixPSlinestyle(pathToFile)

% open loop plots for the benchmark bicycle
figure()
openLoops = bikeData.Medium.openLoops;
num = openLoops.Phi.num;
den = openLoops.Phi.den;
hold all
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
opts.YLim = {[-80, 20], [-540, -90]};
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
closeLeg = legend(lines(8:-1:6), ...
                  {'$\phi$', '$\psi$','$y$'}, ...
                  'Location', 'Southeast', ...
                  'Interpreter', 'Latex');

set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])
filename = 'benchmarkOpen.eps';
pathToFile = ['plots' filesep filename];
saveas(gcf, pathToFile)
fixPSlinestyle(pathToFile)

% handling qualities plot
num = bikeData.Medium.handlingMetric.num;
den = bikeData.Medium.handlingMetric.den;
w = linspace(0.01, 20, 200);
[mag, phase, freq] = bode(tf(num, den), w);
figure()
hold on
metricLine = plot(freq, mag(:)', 'k-', 'Linewidth', 2.0);
ylim([0, 10]);
level1 = line([0, 20], [5, 5]);
level2 = line([0, 20], [8, 8]);
set(level1, 'Color', 'k', 'Linestyle', '--', 'Linewidth', 2.0)
set(level2, 'Color', 'k', 'Linestyle', '--', 'Linewidth', 2.0)
ylabel('Handling Quality Metric')
xlabel('Frequency (rad/sec)')
text(3, 3, 'Level 1')
text(3, 6.5, 'Level 2')
text(3, 9, 'Level 3')
box on
set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])
filename = 'benchmarkHandling.eps';
pathToFile = ['plots' filesep filename];
saveas(gcf, pathToFile)
fixPSlinestyle(pathToFile)

function open_loop_all_bikes(data, linestyles, colors)
% Creates open loop Bode plots of all the bikes.

bikes = fieldnames(data)
freq = {0.1, 20.0};
figure()
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
% used for figure width to height ratio
goldenRatio = (1 + sqrt(5)) / 2;
figWidth = 3.0;
set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])
filename = 'openBode.eps';
pathToFile = ['plots' filesep filename];
print(pathToFile, '-depsc')
fixPSlinestyle(pathToFile)

function handling_all_bikes(data, linestyles, colors)
% Creates handling quality metric for all bikes.
bikes = fieldnames(data);
figure()

% used for figure width to height ratio
goldenRatio = (1 + sqrt(5)) / 2;
figWidth = 3.0;
set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])
w = linspace(0.01, 20, 200);
hold all
speedNames = fieldnames(data.Browser);
fillColors = {[0.98, 0.98, 0.98], [0.93, 0.93, 0.93], [0.85, 0.85, 0.85]};
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

function path_plots(data, linestyles, colors)
% Creates a plot of the path tracking for all bikes at all speeds.

bikes = fieldnames(data);
speedNames = fieldnames(data.Browser);

figure()
goldenRatio = (1 + sqrt(5)) / 2;
figWidth = 3.0;
set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])

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

function plot_io(variable, io, data, linestyles, colors)
% Creates a plot of the time histories of a particular output or input variable
% for three different speeds.
%
% Parameters
% ----------
% variable : string
%   The name of the variable you'd like to plot.
% io : string
%   'input' for input and 'output' for output.
% data : structure
%   Data for a set of bicycles, the first being the benchmark bicycle.
% linestyles : cell array
%   An array of linestyle types, one for each bicycle.
% colors : cell array
%   An array of colors, one for each bicycle.

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
goldenRatio = (1 + sqrt(5)) / 2;
figWidth = 3.0;
set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])


for j = 1:length(speedNames)
    subplot(3, 1, j)
    hold all
    for i = 2:length(bikes)
        oneSpeed = data.(bikes{i}).(speedNames{j});
        time = oneSpeed.time;
        speed = oneSpeed.speed;
        distance = time * speed;
        history = oneSpeed.([io 's'])(:, index);
        plot(distance, history, ...
             'Linestyle', linestyles{i - 1}, ...
             'Color', colors{i - 1})
    end
    xlabel('Distance (m)')
    first = [prettyNames{index} ' ' units{index}];
    second = sprintf(' at %1.1f m/s', speed);
    ylabel([first second], 'Interpreter', 'Latex')
    box on
    xlim([30 190])
    hold off
end
plotAxes = findobj(gcf, 'type', 'axes');
legend(plotAxes(3), bikes(2:end))
filename = [variable '.eps'];
print(['plots' filesep filename], '-depsc')
fixPSlinestyle(['plots' filesep filename])

function phase_portraits(bikeData)
% Creates four phase portrait plots.

figure()
goldenRatio = (1 + sqrt(5)) / 2;
figWidth = 3.0;
set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])

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
