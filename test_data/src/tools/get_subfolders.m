function [subfolders, n_subfolders, sorted_idx] = get_subfolders(directory)
%GET_SUBFOLDERS The function return the subfolders in directory, number of
%subfolders, and sorted index of subfolder created by time.
%
%	See also DIR.

% get subfolders under directory
d = dir(directory);
d = d([d(:).isdir]);

% get rid of '.' and '..'
d = d(~ismember({d(:).name},{'.','..'}));

% get subfolders' name
subfolders = {d(:).name};

% get number of subfolders
n_subfolders = length(d);

% get the latest modified subfolder
[~, Ind] = sort([d(:).datenum], 'descend');

if ~isempty(Ind)
	sorted_idx = Ind;
else
	sorted_idx = 0;
end

end
