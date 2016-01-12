clear; close all; clc;

linestyles = {'-', '-', '-.', ...
              '--', '-.', '--'};
colors = {'k', ...
          [0.5, 0.5, 0.5], ...
          [0.5, 0.5, 0.5], ...
          'k', ...
          'k', ...
          [0.5, 0.5, 0.5]};

trails = linspace(0.1, -0.1, 20);
bikes = {'Browserins', 'Browser', 'Pista', ...
         'Fisher', 'Yellow', 'Yellowrev'};
speed = 3.0;
freqs = linspace(0.01, 20, 100);

% get the directory which this m-file is in
S = dbstack('-completenames');
[CURRENT_DIRECTORY, ~, ~] = fileparts(S(1).file);

for i = 1:length(bikes)
    par = par_text_to_struct([CURRENT_DIRECTORY filesep 'parameters' ...
                              filesep bikes{i} 'Par.txt']);
    for j = 1:length(trails)
        par.c = trails(j);
        [A, B, C, D] = whipple_pull_force_abcd(par, speed);
        % TODO : Don't save the found gains to file.
        data = generate_data(bikes{i}, speed, ...
                             'simulate', false, ...
                             'loopTransfer', false, ...
                             'forceTransfer', {}, ...
                             'fullSystem', false, ...
                             'stateSpace', {A, B, C, D});
        num = data.handlingMetric.num;
        den = data.handlingMetric.den;
        [mag, ~, ~] = bode(tf(num, den), freqs);
        peaks(j) = max(mag);
    end
    plot(trails, peaks, 'Linestyle', linestyles{i}, 'Color', colors{i})
end
plot(trails, ones(size(speeds)) * 5, 'k');
plot(trails, ones(size(speeds)) * 8, 'k');
xlabel('Trail [m]')
ylabel('max(HQW)')
hold off;
