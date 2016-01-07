speeds = linspace(2, 10, 20);
freqs = linspace(0.01, 20, 100);

bikes = {'Benchmark', 'Browserins', 'Browser', 'Pista', ...
         'Fisher', 'Yellow', 'Yellowrev'};
hold on;
for i = 1:length(bikes)
    peaks = zeros(length(speeds));
    for j = 1:length(speeds)
        data = generate_data(bikes{i}, speeds(j));
        num = data.handlingMetric.num;
        den = data.handlingMetric.den;
        [mag, ~, ~] = bode(tf(num, den), freqs);
        peaks(j) = max(mag);
    end
    plot(speeds, peaks)
end
plot(speeds, ones(size(speeds)) * 5);
plot(speeds, ones(size(speeds)) * 8);
hold off;
