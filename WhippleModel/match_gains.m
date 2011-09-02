function bestPar = match_gains(bicycle, speed, stateSpace, tfOutput, experimentTF, guess)
% bestPar = match_gains(bicycle, speed, stateSpace, tfOutput, experimentTF, guess)
%
% Searches for the best set of gains to fit the provided lateral pull force transfer function.
%
% Parameters
% ----------
% bicycle : char
%   The name of the bicycle.
% speed : double
%   The speed of the bicycle.
% tfOutput : char
%   The lateral force transfer function output.
% experimentTF : transfer function
%   The transfer function derived from experimental data.
% guess : matrix, size(1, 6)
%   The guess should contain a starting guess for the five gains and the
%   neuromuscular frequency. [kDelta, kPhiDot, kPhi, kPsi, kY, wnm]
%
% Returns
% -------
% bestPar : matrix, size(1, 6)
%   The parameters that give the best fit of the model to the experimental
%   transfer function.

% find the set of gains and neuro frequency that best fit the transfer
% function
display('hello')
[bestPar, fval, exitflag, output] = ...
    fminsearch(@(x) ...
        transfer_function_mismatch(x, bicycle, speed, experimentTF, stateSpace, tfOutput), ...
        guess, optimset('Display', 'on', 'MaxFunEvals', 20000, 'MaxIter', 10000))
pause

display(sprintf('Function evaluated %d times with %d iterations.', ...
                output.funcCount, output.iterations))

function mismatch = transfer_function_mismatch(parameters, bicycle, ...
    speed, experimentTF, stateSpace, tfOutput)
%
% mismatch = transfer_function_mismatch(parameters, bicycle, speed, ...
%   experimentTF, stateSpace, output)

data = generate_data(bicycle, speed, ...
                     'gains', parameters(1:5), ...
                     'neuroFreq', parameters(6), ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0, ...
                     'stateSpace', stateSpace, ...
                     'forceTransfer', {tfOutput});

modelTF = tf(data.forceTF.(tfOutput).num, ...
             data.forceTF.(tfOutput).den);

% frequency from 1 to 20 radians per second
w = logspace(0, 1.3, 200);

% evaluate the transfer functions at the frequencies
modelResponse = freqresp(modelTF, w);
experimentResponse = freqresp(experimentTF, w);

% calculate the norm of the difference in the two sets of frequencies
mismatch = norm(modelResponse(:) - experimentResponse(:));
% divide by the squareroot of n for the rms
display(sprintf('This is the RMS %1.4f', mismatch/sqrt(200)))

% uncomment to show a graph of the fit
%figure(1)
%clf()
%bode(experimentTF, w)
%hold all
%bode(modelTF, w)
%hold off
