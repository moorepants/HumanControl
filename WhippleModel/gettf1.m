function [g, gd, regmax, err] = gettf1(u, y, nn, tt, flag)
% [g, gd, regmax, err] = gettf1(u, y, nn, tt, flag)
% Returns a transfer function from an input/output time sequence pair and a
% regressor.
%
% Parameters
% ----------
% u : matrix, size(1, n)
%   Input time series.
% y : matrix, size(1, n)
%   Output time series.
% nn : matrix, size(1, 3)
%   Regressor [na nb nk]
% tt : matrix, size(1, n)
%   Time.
% flag : boolean
%   0 to just use the given regressor, and 1 to search for the optimum
%   regressor using the given one as the upper bound.
%
% Returns
% -------
% g : transfer function
%   Estimated continuous transfer function Output/Input (y/u)
% gd : transfer function
%   Estimated discrete transfer function Output/Input (y/u)
% regmax : matrix, size(1, 3)
%   Best regressor when used with optimization (i.e. flag=1)
% err :
%
% Author: Y. Zeyada 16 Aug 2000

% The next 7 lines were added by Ron Hess 11/9/01
ts = tt(6) - tt(5);
[bin,ain] = butter(4, 6 * ts);
[bout,aout] = butter(4, 6 * ts);
in = filter(bin, ain, u);
out = filter(bout, aout, y);
u = in;
y = out;

u = u';
y = y';

NN = nn;

% search for the optimal regressor
if flag==1,
    val = [];
    regr = [];
    % for each order in the denominator of the transfer function
    for NA = 1:nn(1)
        % for each order in the numerator of the transfer function
        for NB = 1:nn(2)
            for NK = -3:nn(3)
                na=NA;
                nb=NB;
                nk=NK;
                % build a * x = b
                j=0;
                a=[];
                b=[];
                for i = 10 + (nk - 1):length(u) - (na + (nk - 1)) - 10,
                    j = j + 1;
                    a(j, 1:na) = y(i:i + na - 1);
                    a(j, na + 1:na + nb) = u(i + na - nk - nb + 1:i + na - nk);
                    b(j, 1) = y(i + na);
                end
                % solve
                x = a \ b;
                xy = x(na:-1:1);
                xu = x(na + nb:-1:na + 1);
                xy1 = zeros(1, nk-1);

                gdd = tf([xu'], [1 -xy' xy1], ts);
                % convert the discrete transfer function to a continous one
                % with a Tustin approximation
                gg = d2c(gdd, 'tustin');
                % cancel any poles and zeros that are close
                gg = minreal(gg);

                % simulate the system with the input and the computed transfer
                % function
                [yc, t, x] = lsim(gg, u, tt);

                % find the error vector between the simulation output and
                % the provided output
                err = yc' - y;
                ems = sqrt(err * err');
                val = [val; 1 / (ems + eps)];
                regr = [regr; NA NB NK];
            end
        end
    end
    maxval = max(val);
    imaxval = find(val == maxval);
    NN = regr(imaxval, :);
end

a = [];
b = [];

if isempty(NN) == 0
    % use the provided regressor
    na = NN(1);
    nb = NN(2);
    nk = NN(3);
    j = 0;
    for i = (10 + (nk - 1)):(length(u) - (na + (nk - 1)) - 10)
        j = j + 1;
        a(j, 1:na) = y(i:(i + na - 1));
        a(j, na + 1:na + nb) = u(i + na-nk-nb + 1:i + na-nk);
        b(j, 1) = y(i + na);
    end
    x = a \ b;
    xy = x(na:-1:1);
    xu = x(na + nb:-1:na + 1);
    xy1 = zeros(1, nk-1);

    gdd = tf([xu'], [1 -xy' xy1 ], ts);
    % find the continuous transfer function using a tustin approximation
    gg = d2c(gdd, 'tustin');
    % cancel any poles and zeros that are close
    gg = minreal(gg);
    % simulate the estimated transfer function with the provided input
    [yc, t, x] = lsim(gg, u, tt);
    %plot(tt, y, tt, yc', 'r')
    %pause;
    err = yc' - y;
    %plot(tt, err)
    % root mean square of y
    rmsy = sqrt(y * y' * ts / max(tt));
    % root mean square error
    rmserr = sqrt(err * err' * ts / max(tt)) / rmsy;
    ems = sqrt(err * err');
    val = 1 / (ems + eps);
    gd = gdd;
    g = gg;
    regmax = NN;
else
    g = [];
    gd = [];
    regmax = [];
end
