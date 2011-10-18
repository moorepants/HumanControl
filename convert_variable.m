function [name, order] = convert_variable(variable, output)
% Returns the name and order of the given variable in the output type.
%
% Parameters
% ----------
% variable : string
%   A variable name.
% output : string.
%   Either `moore`, `meijaard`, `data`.
%
% Returns
% -------
% name : string
%   The variable name in the given output type.
% order : double
%   The order of the variable in the list.

[coordinates, speeds, inputs] = get_variables();

if strcmp(variable(1), 'q')
    type = 'coordinates';
    input = 'moore';
elseif strcmp(variable(1), 'u')
    type = 'speeds';
    input = 'moore';
elseif strcmp(variable(1), 'T')
    type = 'inputs';
    input = 'moore';
elseif regexp(variable, 'Dot$')
    type = 'speeds';
    input = 'meijaard';
elseif regexp(variable, 'Angle$')
    type = 'coordinates';
    input = 'data';
elseif regexp(variable, 'Contact$')
    type = 'coordinates';
    input = 'data';
elseif regexp(variable, 'Rate$')
    type = 'speeds';
    input = 'data';
elseif regexp(variable, 'Torque$')
    type = 'inputs';
    input = 'data';
elseif regexp(variable, 'Force$')
    type = 'inputs';
    input = 'data';
else % search in meijaard coordinates and inputs
    if find(ismember(coordinates.meijaard, variable))
        type = 'coordinates';
        input = 'meijaard';
    elseif find(ismember(inputs.meijaard, variable))
        type = 'inputs';
        input = 'meijaard';
    else
        error('Beep: Done typed yo variable name wrong')
    end
end

if strcmp(type, 'coordinates')
    order = find(ismember(coordinates.(input), variable));
    name = coordinates.(output){order};
elseif strcmp(type, 'speeds')
    order = find(ismember(speeds.(input), variable));
    name = speeds.(output){order};
elseif strcmp(type, 'inputs')
    order = find(ismember(inputs.(input), variable));
    name = inputs.(output){order};
end

function [coordinates, speeds, inputs] = get_variables()

coordinates.moore = {'q1',
                     'q2',
                     'q3',
                     'q4',
                     'q5',
                     'q6',
                     'q7',
                     'q8',
                     'q9',
                     'q10'};

speeds.moore = {'u1',
                'u2',
                'u3',
                'u4',
                'u5',
                'u6',
                'u7',
                'u8',
                'u9',
                'u10'};

inputs.moore = {'T4',
                'T6',
                'T7',
                'F'};

coordinates.data = {'LongitudinalRearContact',
                    'LateralRearContact',
                    'YawAngle',
                    'RollAngle',
                    'PitchAngle',
                    'RearWheelAngle',
                    'SteerAngle',
                    'FrontWheelAngle',
                    'LongitudinalFrontContact',
                    'LateralFrontContact'};

speeds.data = {'LongitudinalRearContactRate',
               'LateralRearContactRate',
               'YawRate',
               'RollRate',
               'PitchRate',
               'RearWheelRate',
               'SteerRate',
               'FrontWheelRate',
               'LongitudinalFrontContactRate',
               'LateralFrontContactRate'};

inputs.data = {'RollTorque',
               'RearWheelTorque',
               'SteerTorque',
               'PullForce'};

coordinates.meijaard = {'xP',
                        'yP',
                        'psi',
                        'phi',
                        'theta',
                        'thetaR',
                        'delta',
                        'thetaF',
                        'xQ',
                        'yQ'};

speeds.meijaard = {'xPDot',
                   'yPDot',
                   'psiDot',
                   'phiDot',
                   'thetaDot',
                   'thetaRDot',
                   'deltaDot',
                   'thetaFDot',
                   'xQDot',
                   'yQDot'};

inputs.meijaard = {'tPhi',
                   'tThetaR',
                   'tDelta',
                   'fB'};
