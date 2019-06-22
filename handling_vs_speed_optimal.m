% Generates a plot of peak HQM versus speed for the optimal bicycles
% presented in the BMD 2016 paper (Figure 10).

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
freqs = linspace(0.01, 40, 200);

speed_strs = {'02', '04', '06', '08', '10'};
bike_lines = zeros(length(speed_strs), 1);

hold on;
for i = 1:length(speed_strs)
    peaks = zeros(1, length(speeds));
    for j = 1:length(speeds)
        data = generate_data(['Optimal', speed_strs{i}], speeds(j), ...
                              'simulate', false, ...
                              'forceTransfer', {}, ...
                              'fullSystem', false);

        lateral_dev_loop = minreal(tf(data.closedLoops.Y.num, data.closedLoops.Y.den));

        if isstable(lateral_dev_loop)
            num = data.handlingMetric.num;
            den = data.handlingMetric.den;
            [mag, ~, ~] = bode(tf(num, den), freqs);
            peaks(j) = max(squeeze(mag));
        else
            peaks(j) = nan;
        end
    end
    bike_lines(i) = plot(speeds, peaks, 'Linestyle', linestyles{i}, 'Color', colors{i});
end
plot(speeds, ones(size(speeds)) * 5, 'k');
plot(speeds, ones(size(speeds)) * 8, 'k');
box on;
xlabel('Speed [m/s]')
ylabel('max(HQM)')
legend(bike_lines, {'2 m/s', '4 m/s', '6 m/s', '8 m/s', '10 m/s'})
ylim([0, 12])
hold off;
