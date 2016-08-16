par = par_text_to_struct('parameters/BenchmarkPar.txt');

guess = zeros(4, 1);
guess(1) = par.c;
guess(2) = par.w;
guess(3) = par.lam;
guess(4) = par.rF;

opts.LBounds = [-inf; 0.010; -pi/2; 1e-10];
opts.UBounds = [inf; inf; pi/2; inf];

[optimal_par, hqm] = cmaes('objective', guess, sqrt(var(guess')'), opts);
