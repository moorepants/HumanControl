function par = par_text_to_struct(pathToFile)
% Returns a structure of the parameters that were stored in a csv text file.
%
% Parameters
% ----------
% pathToFile : string
%   Path to a text file containing the benchmark parameters for a single
%   bicycle. The parameters should be on seperate lines and comma seperated
%   (i.e. c,0.08 or lambda,pi/10)
%
% Returns
% -------
% par : structure
%   A structure containing the bicycle parameters.

fid = fopen(pathToFile);
data = textscan(fid, '%s %s', 'delimiter', ',');
fclose(fid);
names = data{1};
vals = data{2};
for i = 1:length(names)
    par.(names{i}) = str2num(vals{i});
end
