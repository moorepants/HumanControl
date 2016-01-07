speeds = linspace(2, 10, 20);
freqs = linspace(0.01, 20, 100);

bicycle = 'Benchmark'

peaks = zeros(length(speeds));

for i = 1:length(speeds)
    data = generate_data(bicycle, speeds(i));
    num = data.handlingMetric.num;
    den = data.handlingMetric.den;
    [mag, ~, ~] = bode(tf(num, den), freqs);
    peaks[i] = max(mag);
end

plot(speeds, peaks)
