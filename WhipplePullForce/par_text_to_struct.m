function par = par_text_to_struct(pathToFile)
% Returns a structure with the the parameters in a text file.
%
% Parameters
% ----------
% pathToFile : string
%   Path to a text file containing the benchmark parameters for a single
%   bicycle.
% Returns
% -------
% par : structure
%   A structure containing the bicycle parameters.

fid = fopen(pathToFile)
data = textscan(fid, '%s %s', 'delimiter', ',')
fclose(fid)
names = data{1}
vals = data{2}
for i = 1:length(names)
    par.(names{i}) = eval(vals{i});
end
