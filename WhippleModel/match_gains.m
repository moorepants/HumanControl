function modelTF = match_gains(experimentTF)

% get the initial guess for the gains using Ron's technique
data = generate_data('Rigid', 3.2, ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0, ...
                     'forceRollRate', 0);

guess = [data.modelPar.kDelta,
         data.modelPar.kPhiDot,
         data.modelPar.kPhi,
         data.modelPar.kPsi,
         data.modelPar.kY,
         sqrt(data.modelPar.neuroNum)];

% find the set of gains and neuro frequency that best fit the transfer
% function found by the data
bestPar = fminsearch(@(x) transfer_function_mismatch(x, experimentTF), guess);

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

function mismatch = transfer_function_mismatch(parameters, experimentTF)

data = generate_data('Rigid', 3.2, ...
                     'gains', parameters(1:5), ...
                     'neuroFreq', parameters(6), ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0);

modelTF = tf(data.forceRollRateTF.num, data.forceRollRateTF.den);

w = logspace(0, 1.3, 200);

modelResponse = freqresp(modelTF, w);
experimentResponse = freqresp(experimentTF, w);

mismatch = norm(modelResponse(:) - experimentResponse(:));
