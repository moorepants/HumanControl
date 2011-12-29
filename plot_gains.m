function plot_gains(pathToGainFile)
% Plots the gains versus speed.
%
% Parameters
% ----------
% pathToGainFile : string
%   The path to the gain file.

data = importdata(pathToGainFile);

[numSpeeds, numCols] = size(data.data);

figure()
subplot(2, 1, 1)
hold all
plot(data.data(:, 1), data.data(:, 2), '-')
plot(data.data(:, 1), data.data(:, 4), '-')
hold off
ylabel('Gain')
legend('k_\delta', 'k_\phi')

subplot(2, 1, 2)
hold all
plot(data.data(:, 1), data.data(:, 3), '-')
plot(data.data(:, 1), data.data(:, 5), '-')
plot(data.data(:, 1), data.data(:, 6), '-')
hold off
legend('$\dot{\phi}$', 'k_\psi', 'k_{y_d}', 'Interpreter', 'latex')
ylabel('Gain')
xlabel('Speed [m/s]')
