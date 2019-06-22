for speed = [2.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    speed_str = num2str(speed, '%02d');
    path = ['~/Conferences/BMD2016/data/bicycles/Optimal', speed_str, '/Parameters/Optimal', speed_str,  'Benchmark.txt'];
    par = par_text_to_struct(path);
    [A, B, C, D] = whipple_pull_force_abcd(par, speed);
    data = generate_data('Benchmark', speed, 'stateSpace', {A, B, C, D});
    wl = linspace(0.01, 40, 100);
    [mag, phase, freq] = bode(tf(data.handlingMetric.num, data.handlingMetric.den), wl);
    csvwrite(['~/Conferences/BMD2016/data/hqm', speed_str, '.csv'], [wl; squeeze(mag)']')
end

speed_str = num2str(speed, '%02d');
path = '~/Conferences/BMD2016/data/bicycles/Falkorheadonbenchmark/Parameters/FalkorheadonbenchmarkBenchmark.txt';
par = par_text_to_struct(path);
[A, B, C, D] = whipple_pull_force_abcd(par, speed);
data = generate_data('Benchmark', speed, 'stateSpace', {A, B, C, D});
wl = linspace(0.01, 40, 100);
[mag, phase, freq] = bode(tf(data.handlingMetric.num, data.handlingMetric.den), wl);
csvwrite(['~/Conferences/BMD2016/data/hqm-falkor.csv'], [wl; squeeze(mag)']')
