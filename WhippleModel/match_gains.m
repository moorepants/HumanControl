function modelTF = match_gains(experimentTF)
% function modelTF = match_gains(experimentTF)
%
% Searches for the best set of gains to fit the pull force to roll rate
% transfer function.

% get the initial guess for the gains using Ron's technique and the state
% space model
data = generate_data('Rigid', 3.2, ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0);

guess = [data.modelPar.kDelta,
         data.modelPar.kPhiDot,
         data.modelPar.kPhi,
         data.modelPar.kPsi,
         data.modelPar.kY,
         sqrt(data.modelPar.neuroNum)];

stateSpace = {data.modelPar.A,
              data.modelPar.B,
              data.modelPar.C,
              data.modelPar.D};

%w = logspace(0, 1.3, 200);
%figure(1)
%bode(experimentTF, w)
%hold all
%modelTF = tf(data.forceRollRateTF.num, data.forceRollRateTF.den);
%bode(modelTF, w)
%hold off

% find the set of gains and neuro frequency that best fit the transfer
% function found by the data
[bestPar, fval, exitflag, output] = ...
    fminsearch(@(x) transfer_function_mismatch(x, ...
                                               experimentTF, ...
                                               stateSpace), ...
    guess, optimset('MaxFunEvals', 20000, 'MaxIter', 10000));

display(sprintf('Function evaluated %d times with %d iterations.', ...
                output.funcCount, output.iterations))

% now find the transfer function with the best guess
data = generate_data('Rigid', 3.2, ...
                     'gains', bestPar(1:5), ...
                     'neuroFreq', bestPar(6), ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0);

bode(experimentTF)
hold all
modelTF = tf(data.forceRollRateTF.num, data.forceRollRateTF.den);
bode(modelTF)

function mismatch = transfer_function_mismatch(parameters, experimentTF, stateSpace)

data = generate_data('Rigid', 3.2, ...
                     'gains', parameters(1:5), ...
                     'neuroFreq', parameters(6), ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0, ...
                     'stateSpace', stateSpace);

modelTF = tf(data.forceRollRateTF.num, data.forceRollRateTF.den);

w = logspace(0, 1.3, 200);

%figure(1)
%clf()
%bode(experimentTF, w)
%hold all
%bode(modelTF, w)
%hold off

modelResponse = freqresp(modelTF, w);
experimentResponse = freqresp(experimentTF, w);

mismatch = norm(modelResponse(:) - experimentResponse(:));
