bicycle = 'Pista';
par = par_text_to_struct(['parameters/' bicycle 'Par.txt']);

guess = zeros(4, 1);
guess(1) = par.c;
guess(2) = par.w;
guess(3) = par.lam;
guess(4) = par.IFyy;

% wheelbase should always accomdate the mass center
min_wheelbase = (par.mH * par.xH + par.mB * par.xB) / (par.mH + par.mB);

fprintf('The minimum possible wheelbase is %1.4f\n', min_wheelbase);

opts.LBounds = [-inf;          % c
                min_wheelbase; % w
                -pi/2;         % lam
                3e-5];         % IFyy

opts.UBounds = [inf;  % c
                inf;  % w
                pi/2; % lam
                inf]; % IFyy

sigma = [0.5;
         3.0;
         0.3 * pi;
         0.07];
%sigma = sqrt(var(guess')');

[optimal_par, hqm] = cmaes('objective', guess, sigma, opts);
