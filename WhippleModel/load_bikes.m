function data = load_bikes(bikes)
% Returns the data for a set of bicycles.
%
% Parameters
% ----------
% bikes : cell array
%   A cell array that lists the bicyle short names.
%
% Returns
% -------
% data : structure
%   A structure containg a node for each bicycle.

speeds = [2.5, 5.0, 7.5];
speedNames = {'Slow', 'Medium', 'Fast'};

% load the data for all speeds for all the bikes
for i = 1:length(bikes)
    for j = 1:length(speeds)
        data.(bikes{i}).(speedNames{j}) = ...
            generate_data(bikes{i}, speeds(j), 1.0, 0);
    end
end
