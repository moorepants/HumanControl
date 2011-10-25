gains = [76.3808, -0.0516, 7.2456, 0.2632, 0.0708]
neuro = 30

data = generate_data('Rigid', 7.0, 'gains', gains, 'neuroFreq', neuro, 'loopTransfer', 0, 'handlingQuality', 0, 'simulate', 0);

bicycle.A = data.modelPar.A;
bicycle.B = data.modelPar.B;
bicycle.C = data.modelPar.C;
bicycle.D = data.modelPar.D;

bicycle.x = {'xP', 'yP', 'psi', 'phi', 'theta', 'thetaR', 'delta', ...
             'thetaF', 'phiDot', 'thetaRDot', 'deltaDot'};
bicycle.y = {'xP', 'yP', 'psi', 'phi', 'theta', 'thetaR', 'delta', ...
             'thetaF', 'xPDot', 'yPDot', 'psiDot', 'phiDot', ...
             'thetaDot', 'thetaRDot', 'deltaDot', 'thetaFDot', 'xQ', 'yQ'}
bicycle.u = {'tPhi', 'tDelta', 'fB'};

outputs = {'xP', 'yP', 'psi', 'phi', 'theta', 'thetaR', 'delta', ...
           'thetaF', 'xPDot', 'yPDot', 'psiDot', 'phiDot', ...
           'thetaDot', 'thetaRDot', 'deltaDot', 'thetaFDot', ...
           'xQ', 'yQ', 'tDelta'};

inputs = {'fB', 'yc'};

sys = system_state_space(bicycle, gains, neuro, inputs, outputs);

analytic = ss(sys.A, sys.B, sys.C, sys.D, 'StateName', sys.states, ...
    'OutputName', sys.outputs, 'InputName', sys.inputs);

numeric = ss(data.system.A, data.system.B, data.system.C, data.system.D);

figure()
pzplot(analytic, numeric)

% plot the transfer function roll-rate/lateral force for both
figure()
hold all
% plot my analytic model
[num, den] = ss2tf(sys.A, sys.B, sys.C, sys.D, 1);
mine = tf(num(find(strcmp('phiDot', outputs)), :), den)
bode(tf(num(find(strcmp('phiDot', outputs)), :), den))
% plot the data from the simulink model
bode(tf(data.forceTF.PhiDot.num, data.forceTF.PhiDot.den))
[num, den] = ss2tf(data.system.A, data.system.B, data.system.C, data.system.D, 1);
bode(tf(num(12, :), den))

eig(sys.A)
eig(data.system.A)
