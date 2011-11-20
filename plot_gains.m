function plot_gains(pathToGainFile)
% Plots the gains versus speed.
%
% Parameters
% ----------
% pathToGainFile : string
%   The path to the gain file.

data = importdata(pathToGainFile);

[numSpeeds, numCols] = size(data.data);

for i=2:numCols
    subplot(numCols - 1, 1, i - 1)
    plot(data.data(:, 1), data.data(:, i), '.')
    ylabel(data.colheaders{i})
    xlabel('Speed [m/s]')
end
