function data = load_bikes(bikes, input)
% Returns the data for a set of bicycles at three speeds.
%
% Parameters
% ----------
% bikes : cell array
%   A cell array that lists the bicyle short names.
% input : string
%   'Steer' or 'Roll'
%
% Returns
% -------
% data : structure
%   A structure containg a node for each bicycle and each speed.

speeds = [2.5, 5.0, 7.5];
speedNames = {'Slow', 'Medium', 'Fast'};

% load the data for all speeds for all the bikes
for i = 1:length(bikes)
    for j = 1:length(speeds)
        data.(bikes{i}).(speedNames{j}) = ...
            generate_data(bikes{i}, speeds(j), input, 0);
    end
end
