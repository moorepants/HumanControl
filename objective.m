function [peak_hqm] = objective(unknowns)

unknowns

freqs = linspace(0.01, 20, 100);

par = par_text_to_struct('parameters/BenchmarkPar.txt');

par.c = unknowns(1);
par.w = unknowns(2);
par.lam = unknowns(3);
par.rF = unknowns(4);

[A, B, C, D] = whipple_pull_force_abcd(par, 5.0);

data = generate_data('Benchmark', 5.0, ...
                     'simulate', false, ...
                     'forceTransfer', {}, ...
                     'fullSystem', false, ...
                     'stateSpace', {A, B, C, D}, ...
                     'display', 0);

roll_loop = minreal(tf(data.closedLoops.Phi.num, data.closedLoops.Phi.den));

if ~isstable(roll_loop)
    peak_hqm = 100;
else
    num = data.handlingMetric.num;
    den = data.handlingMetric.den;
    [mag, ~, ~] = bode(tf(num, den), freqs);
    peak_hqm = max(mag);
end

peak_hqm
