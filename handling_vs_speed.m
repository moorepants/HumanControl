% Define some linestyles and colors for each of the six bicycles
linestyles = {'-', '-', '-.', ...
              '--', '-.', '--'};
colors = {'k', ...
          [0.5, 0.5, 0.5], ...
          [0.5, 0.5, 0.5], ...
          'k', ...
          'k', ...
          [0.5, 0.5, 0.5]};

speeds = [linspace(3, 4.8, 14), linspace(5, 10, 6)];
freqs = linspace(0.01, 40, 200);  % 0.01 because some have a spike at zero

bikes = {'Browserins', 'Browser', 'Pista', ...
         'Fisher', 'Yellow', 'Yellowrev'};
bike_labels = {'(1)', '(2)', '(3)', '(4)', '(5)', '(6)'};
bike_lines_1 = zeros(1, length(bikes));
bike_lines_2 = zeros(1, length(bikes));

for i = 1:length(bikes)
    peaks = zeros(1, length(speeds));
    lcis = zeros(1, length(speeds));
    for j = 1:length(speeds)
        data = generate_data(bikes{i}, speeds(j), ...
                             'laneType', 'single', ...
                             'simulate', true, ...
                             'loopTransfer', false, ...
                             'forceTransfer', {}, ...
                             'fullSystem', false);
        % Hess handling quality metric
        num = data.handlingMetric.num;
        den = data.handlingMetric.den;
        [mag, ~, ~] = bode(tf(num, den), freqs);
        peaks(j) = max(squeeze(mag));
        % Sandusky's lane change index
        path = data.path;
        time = data.time;
        steer_torque = data.inputs(:, 2);
        roll_rate = data.outputs(:, 12);
        lcis(j) = (max(steer_torque) - min(steer_torque)) / ...
            (max(roll_rate) - min(roll_rate)) / speeds(j);
    end

    subplot(2, 1, 1)
    hold on;
    bike_lines_1(i) = plot(speeds, peaks, 'Linestyle', linestyles{i}, ...
                           'Color', colors{i});
    hold off;

    subplot(2, 1, 2)
    hold on;
    bike_lines_2(i) = plot(speeds, lcis, 'Linestyle', linestyles{i}, ...
                           'Color', colors{i});
    hold off;
end
subplot(2, 1, 1)
hold on;
plot(speeds, ones(size(speeds)) * 5, 'k');
plot(speeds, ones(size(speeds)) * 8, 'k');
box on;
xlabel('Speed [m/s]')
ylabel('max(HQM)')
legend(bike_lines_1, bike_labels)
hold off;
