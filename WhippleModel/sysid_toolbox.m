pathToData = ['..' filesep '..' filesep 'BicycleDataProcessor' filesep 'exports' ...
    filesep 'mat' filesep];

load([pathToData '00264.mat'], 'RollRate', 'SteerRate', 'RollAngle', ...
    'SteerAngle', 'PullForce', 'SteerTorque', 'YawRate', 'YawAngle', ...
    'LateralRearContact')

% build the iddata structure
u = PullForce;

%y = [YawAngle, RollAngle, SteerAngle, RollRate];

y = [YawAngle, RollAngle, SteerAngle, YawRate, RollRate, SteerRate, ...
    LateralRearContact, SteerTorque];

outputNames = {'Yaw Angle', ...
               'Roll Angle', ...
               'Steer Angle', ...
               'Yaw Rate', ...
               'Roll Rate', ...
               'Steer Rate', ...
               'Lateral Rear Contact', ...
               'Steer Torque'};

outputUnits = {'Radian', ...
               'Radian', ...
               'Radian', ...
               'Radian/Second', ...
               'Radian/Second', ...
               'Radian/Second', ...
               'Meter', ...
               'Newton-Meter'}

sampleSize = 1 / 200;

z = iddata(y, u, sampleSize);

set(z, 'InputName', 'Lateral Force')
set(z, 'InputUnit', 'Newtons')

%set(z, 'OutputName', {'Yaw Angle', 'Roll Angle', 'Steer Angle', 'Roll Rate'})
%set(z, 'OutputUnit', {'Radian', 'Radian', 'Radian', 'Radian/Second'})

set(z, 'OutputName', outputNames)
set(z, 'OutputUnit', outputUnits)

ze = z(2000:4000);

% calculate the frequency response
%gs = spa(ze);
%
%bode(gs, 'sd', 3, 'fill')

% find a general model
m = pem(ze)

bode(m)

compare(ze, m)

% go through each input/output pair and get the best parameters for an
% arx model
%naVals = zeros(size(y, 2), 1);
%nbVals = zeros(size(y, 2), 1);
%nkVals = zeros(size(y, 2), 1);
%
%for i = 1:size(y, 2);
    %zs = iddata(y(:, i), u, 1 / 200);
    %set(zs, 'OutputName', outputNames{i}, 'OutputUnit', outputUnits{i})
    %set(zs, 'InputName', 'Lateral Force', 'InputUnit', 'Newtons')
    %naGuess = 6:10;
    %nbGuess = 6:10;
    %nkGuess = delayest(zs(2000:4000));
    %trials = struc(naGuess, nbGuess, nkGuess);
    %arxPar = selstruc(arxstruc(zs(2000:4000), zs, trials), 0);
    %arxMod = arx(zs, arxPar);
    %bode(arxMod)
    %pause
    %naVals(i) = arxPar(1);
    %nbVals(i) = arxPar(2);
    %nkVals(i) = arxPar(3);
%end

%naVals = 10 * ones(size(y, 2))
%nbVals
%nkVals
%
%arxModel = arx(z, 'na', naVals, 'nb', nbVals, 'nk', nbVals)
