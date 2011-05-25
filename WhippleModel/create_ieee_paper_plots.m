function create_ieee_paper_plots(data)
% Creates all of the figures for the IEEE paper.

%loop_shape_example(data.Benchmark)
%open_loop_all_bikes(data)
handling_all_bikes(data)

function loop_shape_example(bikeData)
% Creates the example loop shaping for the bicycle at medium speed.

% used for figure width to height ratio
goldenRatio = (1 + sqrt(5)) / 2;

figure()
freq = {0.1, 20.0};

% the closed delta loop
num = bikeData.Medium.closedLoops.Delta.num;
den = bikeData.Medium.closedLoops.Delta.den;
hold all
deltaBode = bodeplot(tf(num, den), freq);

% the closed phi dot loop
num = bikeData.Medium.closedLoops.PhiDot.num;
den = bikeData.Medium.closedLoops.PhiDot.den;
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

saveas(gcf, 'plots/benchmarkClosed.eps')

% open loop plots for the benchmark bicycle
figure()
num = bikeData.Medium.openLoops.Phi.num;
den = bikeData.Medium.openLoops.Phi.den;
hold all
bodeplot(tf(num, den), freq);

num = bikeData.Medium.openLoops.Psi.num;
den = bikeData.Medium.openLoops.Psi.den;
bodeplot(tf(num, den), freq);

num = bikeData.Medium.openLoops.Y.num;
den = bikeData.Medium.openLoops.Y.den;
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
saveas(gcf, 'plots/benchmarkOpen.eps')

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
saveas(gcf, 'plots/benchmarkHandling.eps')

function open_loop_all_bikes(data)
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

linestyles = {'-', '--', '-.', ...
              '-', '--', '-.'};
colors = {'k', 'k', 'k', ...
          [0.6, 0.6, 0.6], ...
          [0.6, 0.6, 0.6], ...
          [0.6, 0.6, 0.6]};

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
print plots/openBode.eps -depsc

function handling_all_bikes(data)
% Creates handling quality metric for all bikes.
bikes = fieldnames(data);
figure()

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

linestyles = {'-', '--', '-.', ...
              '-', '--', '-.'};
colors = {'k', 'k', 'k', ...
          [0.5, 0.5, 0.5], ...
          [0.5, 0.5, 0.5], ...
          [0.5, 0.5, 0.5]};

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
% used for figure width to height ratio
goldenRatio = (1 + sqrt(5)) / 2;
figWidth = 3.0;
set(gcf, ...
    'PaperSize', [goldenRatio * figWidth, figWidth])
print plots/handling.eps -depsc
fixPSlinestyle('plots/handling.eps')
