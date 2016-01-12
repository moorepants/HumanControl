clear; close all; clc;

% Define some linestyles and colors for each of the six bicycles
linestyles = {'-', '-', '-.', ...
              '--', '-.', '--'};
colors = {'k', ...
          [0.5, 0.5, 0.5], ...
          [0.5, 0.5, 0.5], ...
          'k', ...
          'k', ...
          [0.5, 0.5, 0.5]};

speeds = linspace(3, 10, 20);
freqs = linspace(0.01, 20, 100);

bikes = {'Browserins', 'Browser', 'Pista', ...
         'Fisher', 'Yellow', 'Yellowrev'};
hold on;
for i = 1:length(bikes)
    peaks = zeros(length(speeds));
    for j = 1:length(speeds)
        data = generate_data(bikes{i}, speeds(j), ...
                             'simulate', false, ...
                             'loopTransfer', false, ...
                             'forceTransfer', {}, ...
                             'fullSystem', false);
        num = data.handlingMetric.num;
        den = data.handlingMetric.den;
        [mag, ~, ~] = bode(tf(num, den), freqs);
        peaks(j) = max(mag);
    end
    plot(speeds, peaks, 'Linestyle', linestyles{i}, 'Color', colors{i})
end
plot(speeds, ones(size(speeds)) * 5, 'k');
plot(speeds, ones(size(speeds)) * 8, 'k');
xlabel('Speed [m/s]')
ylabel('max(HQM)')
hold off;
