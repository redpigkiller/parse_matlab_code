function dir_struct = get_subfolders_r(directory, max_recur)
%GET_SUBFOLDERS_R The function return the subfolders in directory, number of
%subfolders, and sorted index of subfolder created by time recursively.
%
% See also DIR, GET_SUBFOLDERS.

if nargin < 2
    max_recur = 10;
end

recur_idx = 0;
dir_struct = pack_info(directory, recur_idx, max_recur);

end

function info_struct = pack_info(directory, recur_idx, max_recur)

info_struct = struct;

if recur_idx >= max_recur
    info_struct.count = 0;
    info_struct.name = [];
    info_struct.array = [];
    return;
end

[subfolders, n_subfolders, sorted_idx] = get_subfolders(directory);

% add field: count
info_struct = setfield(info_struct, 'count', n_subfolders);

% add field: sorted_subfolders name
% add field: sorted_subfolders array
subfolders_name = {};
subfolders_array = [];

for idx = 1 : n_subfolders
    subfolder_name = subfolders{sorted_idx(idx)};
    sub_info_struct = pack_info(fullfile(directory, subfolder_name), ...
        recur_idx + 1, max_recur);

    subfolders_name{end+1} = subfolder_name;
    subfolders_array = [subfolders_array ; sub_info_struct];

end

info_struct = setfield(info_struct, 'name', subfolders_name);
info_struct = setfield(info_struct, 'array', subfolders_array);

end