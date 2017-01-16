function [peak_hqm] = objective(unknowns)

display('==========================================')

unknowns

freqs = linspace(0.01, 40, 200);

par = par_text_to_struct('parameters/BenchmarkPar.txt');

%par.g = 1.6;  % acceleration due to gravity on the moon
%par.g = 3.71;  % acceleration due to gravity on mars
par.c = unknowns(1);
par.w = unknowns(2);
par.lam = unknowns(3);
par.rF = unknowns(4);

speed = 5.0;

[A, B, C, D] = whipple_pull_force_abcd(par, speed);

data = generate_data('Benchmark', speed, ...
                     'simulate', false, ...
                     'forceTransfer', {}, ...
                     'fullSystem', false, ...
                     'stateSpace', {A, B, C, D}, ...
                     'display', 0);

%roll_loop = minreal(tf(data.closedLoops.Phi.num, data.closedLoops.Phi.den));
lateral_dev_loop = minreal(tf(data.closedLoops.Y.num, data.closedLoops.Y.den));

if ~isstable(lateral_dev_loop)
    peak_hqm = max(10, 100 * max(real(pole(lateral_dev_loop))));
else
    num = data.handlingMetric.num;
    den = data.handlingMetric.den;
    [mag, ~, ~] = bode(tf(num, den), freqs);
    peak_hqm = max(mag);
end

peak_hqm
