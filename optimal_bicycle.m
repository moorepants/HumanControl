par = par_text_to_struct('parameters/BenchmarkPar.txt');

guess = zeros(4, 1);
guess(1) = par.c;
guess(2) = par.w;
guess(3) = par.lam;
guess(4) = par.rF;

% wheelbase should always accomdate the mass center
min_wheelbase = (par.mH * par.xH + par.mB * par.xB) / (par.mH + par.mB);

opts.LBounds = [-inf; min_wheelbase; -pi/2; 1e-10];
opts.UBounds = [inf; inf; pi/2; inf];

sigma = [0.5; 3.0; 0.3 * pi; 0.2];
%sigma = sqrt(var(guess')');

[optimal_par, hqm] = cmaes('objective', guess, sigma, opts);
