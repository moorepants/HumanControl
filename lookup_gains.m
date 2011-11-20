function [kDelta, kPhiDot, kPhi, kPsi, kY] = lookup_gains(pathToGains, speed)
% Returns a guess for the gains based on precomputed gains at various speeds
% using linear interpolation/extrapolation.
%
% Parameters
% ----------
% pathToGains : string
%   Path to a csv file with the gains for one bike at various speeds.
% speed : float
%   Speed at which to guess the gains for.
%
% Returns
% -------
% kDelta : float
%   The gain for the steer angle loop.
% kPhiDot : float
%   The gain for the roll rate loop.
% kPhi : float
%   The gain for the roll loop.
% kPsi : float
%   The gain for the heading loop.
% kY : float
%   The gain for the lateral deviation loop.
%
% Notes
% -----
% If there is only one speed in the file, then the function returns the value
% of those gains no matter which speed you specify.

try
    contents = importdata(pathToGains);
catch err
    display(sprintf('There is no gain file: %s.', pathToGains))
    rethrow(err)
end

speeds = contents.data(:, 1);

if length(speeds) == 1
    guesses = contents.data(2:end);
    display_if(sprintf(['Gain guess may be bad, please provide more than one ' ...
                    'speed of gain guesses in %s'], pathToGains))
else
    guesses = zeros(5);
    % interpolate/extrapolate the data
    for i = 2:6
        gains = contents.data(:, i);
        % use extrapolation if needed
        if speed > speeds(end) || speed < speeds(1)
            guesses(i - 1) = interp1(speeds, gains, speed, 'linear', 'extrap');
        else
            guesses(i - 1) = interp1(speeds, gains, speed, 'linear');
        end
    end
end

kDelta = guesses(1);
kPhiDot = guesses(2);
kPhi = guesses(3);
kPsi = guesses(4);
kY = guesses(5);
