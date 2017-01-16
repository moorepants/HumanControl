par = 'w';
start = 0.1;
stop = 1.1;
num = 30;

% Generate some parameter values to compute the HQM at.
par_vals = linspace(start, stop, num);
peaks = zeros(1, length(par_vals));
% These are frequencies that the HQM transfer function should be evaluated
% at.
freqs = linspace(0.01, 20, 100);

for j = 1:length(par_vals)
    % Load the benchmark bicycle parameters from the file.
    par = par_text_to_struct('parameters/BenchmarkPar.txt');
    % Change the parameter you are plotting to the jth value.
    par.(par) = par_vals(j);
    % Generate the state space form of the plant (bicycle) at 5.0 m/s and
    % the specified physical parameters.
    [A, B, C, D] = whipple_pull_force_abcd(par, 5.0);
    % Generate the HQM transfer function with the specified state space
    % matrices.
    data = generate_data('Benchmark', 5.0, ...
                         'simulate', false, ...
                         'loopTransfer', false, ...
                         'forceTransfer', {}, ...
                         'fullSystem', false, ...
                         'stateSpace', {A, B, C, D});
    % Find the magnitude of the HQM transfer function at all of the desired
    % freqencies.
    num = data.handlingMetric.num;
    den = data.handlingMetric.den;
    [mag, ~, ~] = bode(tf(num, den), freqs);
    % Store the peak value of the HQM.
    peaks(j) = max(mag);
end

% Plot the parameter values versus the peak HQM.
plot(w, peaks)
xlabel(par)
ylabel('Max of HQM')
